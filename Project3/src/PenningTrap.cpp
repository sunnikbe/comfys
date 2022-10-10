#include "PenningTrap.hpp"
#include "Particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0d_in)
{
  B0_ = B0_in;
  V0d_ = V0d_in;
}

// Methods returning member variable values
// Add a new particle to the Penning trap
void PenningTrap::add_particle(Particle particle_in)
{
  particles_.push_back(particle_in);
}

// External electric field at point r
arma::vec PenningTrap::external_E_field(arma::vec r)
{
  arma::vec tmp = {1,1,-2};
  return V0d_*r*tmp;
}

// External magneticfield at point r
arma::vec PenningTrap::external_B_field(arma::vec v)
{
  arma::vec B = {0,0,B0_};
  return arma::cross(v,B);
}

// force on particle i from particle j
arma::vec PenningTrap::force_particle(int i, int j)
{
  double k_e = 1.38935333*std::pow(10,5);
  double q_j = particles_.at(j).charge();
  arma::vec r_i = particles_.at(i).position();
  arma::vec r_j = particles_.at(j).position();
  return k_e*q_j*(r_i-r_j)/std::pow(norm(r_i-r_j),3);
}

// total force on particle i from external forces
arma::vec PenningTrap::total_force_external(int i)
{
  double q_i = particles_.at(i).charge();
  arma::vec r_i = particles_.at(i).position();
  arma::vec v_i = particles_.at(i).velocity();
  return q_i*external_E_field(r_i) + q_i*external_B_field(v_i);
}

// total force on particle i from other particles
arma::vec PenningTrap::total_force_particle(int i)
{
  arma::vec F_particle = {0,0,0};
  // need more work
  for (int j = 0; j < particles_.size(); j++)
  {
    if (j != i)
    {
      F_particle += force_particle(i,j);
    }
  }
  return F_particle;
}

// Total force on particle i
arma::vec PenningTrap::total_force(int i)
{
  return total_force_external(i) + total_force_particle(i);
}

// Evolve PenningTrap in time (Problem 7)
// Runge-Kutta:


// Forward Euler for particle i in PenningTrap:
void PenningTrap::evolve_fEuler(int i, arma::vec t, double h)
{
  // From particles when put into penning trap
  arma::vec r_i_0 = particles_.at(i).position();
  arma::vec v_i_0 = particles_.at(i).velocity();
  double q_i = particles_.at(i).charge();
  double m_i = particles_.at(i).mass();

  // for f function
  double x_0 = r_i_0(0);
  double z_0 = r_i_0(2);
  double v_0 = v_i_0(1);
  double w_0 = (q_i*B0_)/m_i;
  double w_z_squared = ((2.*q_i)/(m_i))*V0d_;
  double w_z = sqrt(w_z_squared);

  double w_p = w_0/2. + sqrt(w_0*w_0 - 2*w_z_squared);
  double w_m = w_0/2. - sqrt(w_0*w_0 - 2*w_z_squared);

  double A_p = (v_0 + w_m*x_0)/(w_m - w_p);
  double A_m = -(v_0 + w_p*x_0)/(w_m - w_p);

  std::complex<double> J = sqrt(-1);

  // Forward Euler:
  arma::vec r_i = arma::vec(t.size()).zeros();
  arma::vec v_i = arma::vec(t.size()).zeros();

  // Initial conditions
  r_i.col(0) = r_i_0;
  v_i.col(0) = v_i_0;

  // Forward Euler algorithm
  for (int k = 1; k < t.size(); k++)
  {
    // f function
    std::complex<double> f = A_p*exp(-J*w_p*t(k)) + A_m*exp(-J*w_m*t(k));

    std::complex<double> r_x = f.real();
    std::complex<double> r_y = f.imag();
    double r_z = z_0*cos(w_z*t(k));

    // Y function, position vector
    arma::vec Y = arma::vec("r_x r_y r_z");

    v_i.col(k + 1) = v_i.col(k) + h*Y;
    r_i.col(k + 1) = r_i.col(k) + h*v_i.col(k);
  }

  printf("Velocity:\n");
  printf("Position:\n");
  v_i.print(std::cout);
  r_i.print(std::cout);
}
