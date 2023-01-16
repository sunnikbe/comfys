// The Particle class (Problem 5)

#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>
#include <string>
#include <cmath>


class Particle
{
  public:
    double q_, m_; // charge, mass
    arma::vec r_, v_; // position, velocity

    // Constructor
    Particle(double, double, arma::vec, arma::vec);

    // Methods
    double charge();
    double mass();
    arma::vec position();
    arma::vec velocity();

    // Method giving info
    std::string info();
};

#endif
