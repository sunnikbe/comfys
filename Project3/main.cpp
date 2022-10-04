#include "PenningTrap.hpp"

int main()
{
  // Make a particle of type Particle(class) to test
  Particle particle_1 = Particle(1.0, 0.2, arma::vec("1 2 3"), arma::vec("3 2 1"));
  // printing info to test
  std::cout << particle_1.info();

  PenningTrap trap = PenningTrap(96.5, 9.65*std::pow(10,8), std::pow(10,4));
  trap.add_particle(particle_1);

  Particle particle_2 = Particle(1.0, 0.1, arma::vec("3 2 1"), arma::vec("1 2 3"));
  trap.add_particle(particle_2);

  arma::vec r_test,ext_force,part_force;
  //r_test = trap.external_B_field(particle_1.position());
  r_test = trap.total_force(0);
  r_test.print(std::cout);


  return 0;
}
