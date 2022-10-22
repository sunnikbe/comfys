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

// External electric field at point r
vec PenningTrap::external_E_field(double t, vec r, double f, double omega_v)
{
  vec tmp = {1,1,-2};
  vec E_field;
  //std::cout << V0d_*r%tmp << std::endl;
  if (time_dependent_potential_ == 1 && norm(r) > d_)
  {
    E_field = V0d_*(1 + f*std::cos(omega_v*t))*r%tmp;
  }
  else if (time_dependent_potential_ == 0 && norm(r) > d_)
  {
    E_field = V0d_*r%tmp;
  }
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
  if (norm(r) > d_)
  {
    B_field = cross(v,B);
  }
  else{
    B_field = {0,0,0};
  }
  //std::cout << cross(v,B) << std::endl;
  return B_field;
}

// force on particle i from particle j
vec PenningTrap::force_particle(int i, int j, vec r_i)
{
  double k_e = 1.38935333*std::pow(10,5);
  double q_j = particles_.at(j).charge();
  vec r_j = particles_.at(j).position();
  return k_e*q_j*(r_i-r_j)/std::pow(norm(r_i-r_j),3);
}

// total force on particle i from external forces
vec PenningTrap::total_force_external(int i, double t, vec r_i, vec v_i, double f, double omega_v)
{
  double q_i = particles_.at(i).charge();
  return q_i*external_E_field(t,r_i,f,omega_v) + q_i*external_B_field(v_i, r_i);
}

// total force on particle i from other particles
vec PenningTrap::total_force_particle(int i, vec r)
{
  vec F_particle = {0,0,0};
  double q_i = particles_.at(i).charge();
  for (int j = 0; j < particles_.size(); j++)
  {
    if (j != i)
    {
      F_particle += q_i*force_particle(i,j,r);
    }
  }
  return F_particle;
}

// Total force on particle i
vec PenningTrap::total_force(int i, double t, vec r, vec v, double f, double omega_v)
{
  vec F_tot;
  if (use_Coulomb_interactions_ == 1)
  {
    F_tot = total_force_external(i,t,r,v,f,omega_v) + total_force_particle(i,r);
  }
  else
  {
    F_tot = total_force_external(i,t,r,v,f,omega_v);
  }
  return F_tot;
}

// Evolve PenningTrap in time (Problem 7)
// Runge-Kutta:
void PenningTrap::evolve_RK4(double t, double dt, double f, double omega_v)
{
  for (int i = 0; i < particles_.size(); i++)
  {
    double q = particles_.at(i).charge();
    double m = particles_.at(i).mass();
    vec v = particles_.at(i).velocity();
    vec r = particles_.at(i).position();

    vec k1_v = dt*total_force(i,t,r,v,f,omega_v);
    vec k1_r = dt*v;
    
    vec k2_r = dt*(v + 0.5*k1_r);
    vec k2_v = dt*total_force(i,t+0.5*dt,r+0.5*k1_r,v+0.5*k1_v,f,omega_v);

    vec k3_r = dt*(v + 0.5*k2_r);
    vec k3_v = dt*total_force(i,t+0.5*dt,r+0.5*k2_r,v+0.5*k2_v,f,omega_v);

    vec k4_r = dt*(v + k3_r);
    vec k4_v = dt*total_force(i,t+dt,r+k3_r,v+k3_v,f,omega_v);
    
    vec new_r = r + 1./6*(k1_r + 2*k2_r + 2*k3_r + k4_r); 
    vec new_v = v + 1./6*(k1_v + 2*k2_v + 2*k3_v + k4_v)/m;

    Particle updated_ = Particle(q, m ,new_r, new_v);

    particles_.at(i) = updated_;
    //std::cout << "Particle" << i << "\n" << particles_.at(i).info();
  }
}

// Forward Euler for all particles in PenningTrap:
void PenningTrap::evolve_fEuler(double t, double dt, double f, double omega_v)
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
    vec F = total_force(i,t,r,v,f,omega_v);

    vec new_v = particles_.at(i).velocity() + dt*F/m;
    vec new_r = particles_.at(i).position() + dt*particles_.at(i).velocity();

    Particle updated_ = Particle(q, m, new_r, new_v);

    particles_.at(i) = updated_;
    std::cout << "Particle " << i << ":\n" << particles_.at(i).info();
  }
}

// write the positions into file
void PenningTrap::write_to_file(std::string filename)
{
  // append new content to file instead of overwriting it
  std::ofstream ofile(filename, std::ios::out | std::ios::app);
  for (int i = 0; i < particles_.size(); i++)
  {
    vec r = particles_.at(i).position();
    ofile << std::setw(16) << std::scientific << r(0)
    << std::setw(16) << std::scientific << r(1)
    << std::setw(16) << std::scientific << r(2);
  }
  ofile << std::endl;
}
