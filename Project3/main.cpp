#include "PenningTrap.hpp"

int main()
{
  // Make a particle of type Particle(class) to test
  Particle particle_1 = Particle(1.0, 0.1, arma::vec("1 1 1"), arma::vec("1 1 1"));
  // printing info to test
  // std::cout << particle_1.info();

  // create an instance of the PenningTrap class with name "trap"
  PenningTrap trap = PenningTrap(96.5, 9.65);
  // add particle_1 into the trap
  trap.add_particle(particle_1);

  // create a second particle and add it into the trap
  Particle particle_2 = Particle(2.0, 0.2, arma::vec("2 2 2"), arma::vec("2 2 2"));
  trap.add_particle(particle_2);

  // // test values received from the functions in the PenningTrap class
  // // change the function name if needed to test the function
  // arma::vec r_test;
  // // test function that takes position vector inpuy
  // //r_test = trap.external_B_field(particle_1.position());
  //
  // // test function that takes particle index input
  // r_test = trap.total_force(0);
  // r_test.print(std::cout);

  // Testing the forward Euler:
//---------------------------------------
  int n = 4; // number of steps
  int N = n - 1; // interior points
  double dt = 1./n; // stepsize
  arma::vec t = arma::linspace(0, n, N); // makes t values with stepsize dt
  arma::vec timesteps = arma::vec(N).fill(dt); // evolves fEuler N times.

  arma::vec Y = arma::vec("10 3 2");

  for (int i = 0; i <N; i++)
  {
    double f = t(i);
    printf("t = %f:\n", f);
    trap.evolve_fEuler(timesteps(i), Y);
  }
  return 0;
}
