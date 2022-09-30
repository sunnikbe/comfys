#include "PenningTrap.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particles_in);
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  particles_ = particles_in;
}

// Methods returning member variable values
std::vector<Particle> PenningTrap::particles()
{
  return particles_;
}
