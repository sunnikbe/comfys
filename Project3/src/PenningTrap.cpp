#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0d_in, double d_in, 
bool time_dependent_potential_in, bool use_Coulomb_interactions_in)
{
  B0_ = B0_in;
  V0d_ = V0d_in;
  d_ = d_in;
  time_dependent_potential_ = time_dependent_potential_in;
  use_Coulomb_interactions_ = use_Coulomb_interactions_in;
}

// Methods returning member variable values
// Add a new particle to the Penning trap
void PenningTrap::add_particle(Particle particle_in)
{
  particles_.push_back(particle_in);
}

// update a given particle inside the trap with a new Particle class
void PenningTrap::update_particle(Particle particle_in, int i)
{
  particles_.at(i) = particle_in;
}

// function that turns on or off the time dependent potential
// false = off ; true = on
void PenningTrap::time_dependent_potential(bool statement)
{
  time_dependent_potential_ = statement;
}

// function that turnson or off the Coulomb interactions
// false = off ; true = on
void PenningTrap::use_Coulomb_interactions(bool statement)
{
  use_Coulomb_interactions_ = statement;
}

// External electric field at point r
vec PenningTrap::external_E_field(double t, vec r, double f, double omega_v)
{
  // tmp vector to set V(x,y,z) = V0d*(x + y - 2z)
  vec tmp = {1,1,-2};
  vec E_field;
  // check if time dependency is on
  // or if the particle is still in the trap
  if (time_dependent_potential_ == 1 && norm(r) < d_)
  {
    E_field = V0d_*(1 + f*std::cos(omega_v*t))*r%tmp;
  }
  else if (time_dependent_potential_ == 0 && norm(r) < d_)
  {
    E_field = V0d_*r%tmp;
  }
  // if not in the trap, set the electric field to zero
  else
  {
    E_field = {0,0,0};
  }
  return E_field;
}

// External magneticfield at point r
vec PenningTrap::external_B_field(vec v, vec r)
{
  vec B_field;
  vec B = {0,0,B0_};
  // check if the particle is still in the trap
  if (norm(r) < d_)
  {
    B_field = cross(v,B);
  }
  // if not set the magnetic field to zero
  else{
    B_field = {0,0,0};
  }
  return B_field;
}

// force on particle i from particle j
vec PenningTrap::force_particle(int i, int j)
{
  double k_e = 1.38935333*std::pow(10,5);
  double q_j = particles_.at(j).charge();
  vec r_i = particles_.at(i).position();
  vec r_j = particles_.at(j).position();
  return k_e*q_j*(r_i-r_j)/std::pow(norm(r_i-r_j),3);
}

// total force on particle i from external forces
vec PenningTrap::total_force_external(int i, double t, double f, double omega_v)
{
  double q_i = particles_.at(i).charge();
  vec r_i = particles_.at(i).position();
  vec v_i = particles_.at(i).velocity();
  return q_i*external_E_field(t,r_i,f,omega_v) + q_i*external_B_field(v_i,r_i);
}

// total force on particle i from other particles
vec PenningTrap::total_force_particle(int i)
{
  vec F_particle = {0,0,0};
  double q_i = particles_.at(i).charge();
  for (int j = 0; j < particles_.size(); j++)
  {
    if (j != i)
    {
      F_particle += q_i*force_particle(i,j);
    }
  }
  return F_particle;
}

// Total force on particle i
vec PenningTrap::total_force(int i, double t, double f, double omega_v)
{
  vec F_tot;
  if (use_Coulomb_interactions_ == 1)
  {
    F_tot = total_force_external(i,t,f,omega_v) + total_force_particle(i);
  }
  else
  {
    F_tot = total_force_external(i,t,f,omega_v);
  }
  return F_tot;
}

// Evolve PenningTrap in time (Problem 7)
// Runge-Kutta:
void PenningTrap::evolve_RK4(double dt, double t, double f, double omega_v)
{
  int n = particles_.size();
  vec q(n);
  vec m(n);
  mat r(3,n);
  mat v(3,n);
  for (int i = 0; i < n; i++)
  {
    q(i) = particles_.at(i).charge();
    m(i) = particles_.at(i).mass();
    r.col(i) = particles_.at(i).position();
    v.col(i) = particles_.at(i).velocity();
  }

  mat k1_v(3,n);
  mat k1_r(3,n);
  // compute k1_r and k1_v
  for (int i = 0; i < particles_.size(); i++)
  {
    k1_v.col(i) = dt*total_force(i,t,f,omega_v)/m(i);
    k1_r.col(i) = dt*v.col(i);
  }
  // update pos and vel for particles corresponding to k1_r and k1_v
  for (int i = 0; i < n; i++)
  {
    vec v_k1 = particles_.at(i).velocity() + 0.5*k1_v.col(i);
    vec r_k1 = particles_.at(i).position() + 0.5*k1_r.col(i);
    particles_.at(i) = Particle(q(i),m(i),r_k1,v_k1);
  }

  mat k2_v(3,n);
  mat k2_r(3,n);
  // compute k2_r and k2_v
  for (int i = 0; i < particles_.size(); i++)
  {
    k2_v.col(i) = dt*total_force(i,t+0.5*dt,f,omega_v)/m(i);
    k2_r.col(i) = dt*(v.col(i) + 0.5*k1_v.col(i));
  }
  // update pos and vel for particles corresponding to k2_r and k2_v
  for (int i = 0; i < n; i++)
  {
    vec v_k2 = particles_.at(i).velocity() + 0.5*k2_v.col(i);
    vec r_k2 = particles_.at(i).position() + 0.5*k2_r.col(i);
    particles_.at(i) = Particle(q(i),m(i),r_k2,v_k2);
  }

  mat k3_v(3,n);
  mat k3_r(3,n);
  // compute k3_r and k3_v
  for (int i = 0; i < particles_.size(); i++)
  {
    k3_v.col(i) = dt*total_force(i,t+0.5*dt,f,omega_v)/m(i);
    k3_r.col(i) = dt*(v.col(i) + 0.5*k2_v.col(i));
  }
  // update pos and vel for particles corresponding to k3_r and k3_v
  for (int i = 0; i < n; i++)
  {
    vec v_k3 = particles_.at(i).velocity() + k3_v.col(i);
    vec r_k3 = particles_.at(i).position() + k3_r.col(i);
    particles_.at(i) = Particle(q(i),m(i),r_k3,v_k3);
  }

  mat k4_v(3,n);
  mat k4_r(3,n);
  // compute k4_r and k4_v
  for (int i = 0; i < particles_.size(); i++)
  {
    k4_v.col(i) = dt*total_force(i,t+dt,f,omega_v)/m(i);
    k4_r.col(i) = dt*(v.col(i) + k3_v.col(i));
  }

  vec new_r, new_v;
  // update new pos and vel for particles given ki_r and ki_v
  for (int i = 0; i < n; i++)
  {
    new_r = r.col(i) + 1./6*(k1_r.col(i) + 2*k2_r.col(i) + 2*k3_r.col(i) + k4_r.col(i));
    new_v = v.col(i) + 1./6*(k1_v.col(i) + 2*k2_v.col(i) + 2*k3_v.col(i) + k4_v.col(i));
    particles_.at(i) = Particle(q(i),m(i),new_r,new_v);
  }
}

// Forward Euler for all particles in PenningTrap:
void PenningTrap::evolve_fEuler(double dt, double t, double f, double omega_v)
{
  // Evolves v and r for all particles one timestep dt
  // f has to be a function where
  // Y(t) = [Re(f(t)), Im(f(t)), z(t)]

  for (int i = 0; i < particles_.size(); i++)
  {
    double q = particles_.at(i).charge();
    double m = particles_.at(i).mass();
    vec v = particles_.at(i).velocity();
    vec r = particles_.at(i).position();
    vec F = total_force(i,t,f,omega_v);

    vec new_v = particles_.at(i).velocity() + dt*F/m;
    vec new_r = particles_.at(i).position() + dt*particles_.at(i).velocity();

    Particle updated_ = Particle(q, m, new_r, new_v);

    particles_.at(i) = updated_;
    //std::cout << "Particle " << i << ":\n" << particles_.at(i).info();
  }
}

// write the positions and velocities to file
// if only one statement is given it will only write to that file
void PenningTrap::write_to_file(std::string filename_pos, std::string filename_vel)
{
  // append new content to file instead of overwriting it
  if (filename_pos != "none")
  {
    std::ofstream ofile_pos(filename_pos, std::ios::out | std::ios::app);
    for (int i = 0; i < particles_.size(); i++)
    {
      vec r = particles_.at(i).position();
      ofile_pos << std::setw(16) << std::scientific << r(0)
      << std::setw(16) << std::scientific << r(1)
      << std::setw(16) << std::scientific << r(2);
    }
    ofile_pos << std::endl;
  }

  // write the velocities to file
  if (filename_vel != "none")
  {
    std::ofstream ofile_vel(filename_vel, std::ios::out | std::ios::app);
    for (int i = 0; i < particles_.size(); i++)
    {
      vec v = particles_.at(i).velocity();
      ofile_vel << std::setw(16) << std::scientific << v(0)
      << std::setw(16) << std::scientific << v(1)
      << std::setw(16) << std::scientific << v(2);
    }
    ofile_vel << std::endl;
  }
}

// function that checks how many particles are inside the trap
// at the moment of the call
int PenningTrap::particles_in_trap()
{
  int particles_inside = 0;
  for (int i = 0; i < particles_.size(); i++)
  {
    vec r_i = particles_.at(i).position();
    if (norm(r_i) < d_)
    {
      particles_inside += 1;
    }
  }
  return particles_inside;
}
