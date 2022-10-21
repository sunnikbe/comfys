// The PenningTrap class (Problem 6)

#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <Particle.hpp>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>


class PenningTrap
{
public:
  double B0_; // magnetic field strength
  double V0d_; // applied potential / characteristic dimension^2
  double d_;
  bool time_dependent_potential_;
  bool use_Coulomb_interactions_;
  std::vector<Particle> particles_;

  // Constructor
  PenningTrap(double B0_in, double V0d_in, double d_in, 
  bool time_dependent_potential_in, bool use_Coloumb_interactions_in);

  // Methods
  void add_particle(Particle particle_in);

  arma::vec external_E_field(double t, arma::vec r, double f, double omega_v);

  arma::vec external_B_field(arma::vec v, arma::vec r);

  arma::vec force_particle(int i, int j, arma::vec r_i);

  arma::vec total_force_external(int i, double t, arma::vec r, arma::vec v, double f, double omega_v);

  arma::vec total_force_particle(int i, arma::vec r);
  
  arma::vec total_force(int i, double t, arma::vec r, arma::vec v, double f, double omega_v);

  void evolve_RK4(double t, double dt, double f = 0, double omega_v = 0);

  void evolve_fEuler(double t, double dt, double f = 0, double omega_v = 0);

  void write_to_file(std::string filename);

};

#endif
