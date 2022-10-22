#include "PenningTrap.hpp"

int main()
{
  bool time_dependent_potential = false;
  bool use_Coulomb_interactions = true;

  arma::vec r,v;
  double B0 = 96.5; // magnetic field 
  double d = 500; // size of trap
  double V0d = 9.65; // V0/d^2
  double m = 40.078; // mass of particle in atomic mass
  double q = 1.0; // charge of particle in elementary charge
  double domega = 0.02;
  arma::vec freq = {0.1,0.4,0.7};
  arma::vec omega_v = arma::linspace(0.2,2.5,(2.5-0.2)/domega);
  // create an instance of the PenningTrap class with name "trap"
  PenningTrap trap = PenningTrap(B0, V0d, d, time_dependent_potential,
  use_Coulomb_interactions);

  // problem 8 for testing the PenningTrap for two particles
  // adding M particles to the trap
  r = {20,0,20};
  v = {0,25,0};
  Particle particle_1 = Particle(q,m,r,v);
  trap.add_particle(particle_1);

  r = {25,25,0};
  v = {0,40,5};
  Particle particle_2 = Particle(q,m,r,v);
  trap.add_particle(particle_2);

  int T = 50;
  double dt = 1e-3;
  int N = T/dt;
  arma::vec t = arma::linspace(0, T, N); // makes t values with stepsize dt
  arma::vec timesteps = arma::vec(N).fill(dt); // evolves fEuler N times.
  std::string filename = "RK4.txt";
  //std::string filename = "fEuler.txt";
  // clears previous content from the file
  std::ofstream ofile(filename);

  trap.write_to_file(filename);
  for (int i = 1; i <N; i++)
  {
    double f = t(i);
    trap.evolve_RK4(t(i), timesteps(i));
    trap.write_to_file(filename);
  }

  // problem 9, not implemented correctly
  int M = 0;
  for (int i = 0; i < M; i++)
  {
    r = arma::vec(3).randn()*0.1*d;
    v = arma::vec(3).randn()*0.1*d;
    // create a particle of the Particle class
    Particle new_particle = Particle(q, m, r, v);
    //std::cout << new_particle.info();
    // add the created particle into the Penning trap
    trap.add_particle(new_particle);
  }

  //for (int j = 0; j < freq.size(); j++)
  //{
    //for (int k = 0; k < omega_v.size(); k++)
    //{

      // test values received from the functions in the PenningTrap class
      // change the function name if needed to test the function
      //arma::vec r_test;
      //test function that takes position vector inpuy
      //r_test = trap.external_B_field(particle_1.position());
      //
      //test function that takes particle index input
      //r_test = trap.total_force(0);
      //r_test.print(std::cout);

    //}
  //}
  return 0;
}
