// The PenningTrap class (Problem 6)

#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <Particle.hpp>
#include <vector>
#include <cmath>
#include <complex>


class PenningTrap
{
public:
  double B0_; // magnetic field strength
  double V0d_; // applied potential / characteristic dimension^2
  std::vector<Particle> particles_;

  // Constructor
  PenningTrap(double B0_in, double V0d_in);

  // Methods
  void add_particle(Particle particle_in);

  arma::vec external_E_field(arma::vec r);

  arma::vec external_B_field(arma::vec v);

  arma::vec force_particle(int i, int j);

  arma::vec total_force_external(int i);

  arma::vec total_force_particle(int i);

  arma::vec total_force(int i);

  void evolve_fEuler(double dt);

};

#endif
