#include "Particle.hpp"

int main()
{
  // Make a particle of type Particle(class) to test
  Particle particle_1 = Particle(1.0, 0.2, arma::vec("1 2 3"), arma::vec("3 2 1"));
  // printing info to test
  std::cout << particle_1.info();

  return 0;
}
