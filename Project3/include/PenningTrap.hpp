// The PenningTrap class (Problem 6)

#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <Particle.hpp>
#include <vector>


class PenningTrap
{
public:
  double B0_; // magnetic field strength
  double V0_; // applied potential
  double d_; // characteristic dimension
  std::vector<Particle> particles_;

  // Constructor
  PenningTrap(double, double, double, std::vector<Particle>);

  // Methods
  std::vector<Particle> particles();

};

#endif
