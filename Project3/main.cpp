#include "PenningTrap.hpp"
#include <time.h>

int main()
{
  // declarations of vatiables used for the Penningtrap
  arma::vec r,v;
  double B0 = 96.5; // magnetic field 
  double d = 500; // size of trap
  double V0d = 9.65; // V0/d^2
  double m = 40.078; // mass of particle in atomic mass
  double q = 1.0; // charge of particle in elementary charge
  double domega = 0.02;
  // create an instance of the PenningTrap class with name "trap"
  // statements for time dependent potential and Coulomb interactions are optional
  PenningTrap trap = PenningTrap(B0, V0d, d);

  // problem 8 for testing the PenningTrap for one and two particles
  // adding particle 1 to the trap
  r = {20,0,20};
  v = {0,25,0};
  Particle particle_1 = Particle(q,m,r,v);
  trap.add_particle(particle_1);

  // declaring how long we shall evolve the particles inside the trap
  int T = 50;
  double dt = 1e-3;
  int N = T/dt;
  arma::vec t = arma::linspace(0, T, N); // makes t values with stepsize dt

  // filenames created when running the write_to_file function in the PenningTrap class
  std::string filename_pos = "RK4_pos_single_particle.txt";
  //std::string filename = "fEuler.txt";
  std::ofstream ofile_pos;
  std::ofstream ofile_vel;
  // clears previous content from the file
  // this is needed as the write_to_file function uses append to write to file
  ofile_pos.open(filename_pos, std::ofstream::out | std::ofstream::trunc);

  // write down the starting values for position and velocity
  trap.write_to_file(filename_pos);
  for (int i = 1; i < N; i++)
  {
    // evolves the particles with Runge-Kutta4
    trap.evolve_RK4(dt);
    trap.write_to_file(filename_pos);
  }
  ofile_pos.close();

  // after the end of the loop we need to update the particles inside the trap
  // if we want to continue with the same trap instance
  particle_1 = Particle(q,m,r,v);
  trap.update_particle(particle_1,0);

  // evolves the same particle with forward Euler
  filename_pos = "fEuler_pos.txt";
  ofile_pos.open(filename_pos, std::ofstream::out | std::ofstream::trunc);

  trap.write_to_file(filename_pos);
  for (int i = 1; i < N; i++)
  {
    trap.evolve_fEuler(dt);
    trap.write_to_file(filename_pos);
  }
  ofile_pos.close();

  // two particle test

  particle_1 = Particle(q,m,r,v);
  trap.update_particle(particle_1,0);
  
  // adds a second particle to the trap
  r = {25,25,0};
  v = {0,40,5};
  Particle particle_2 = Particle(q,m,r,v);
  trap.add_particle(particle_2);

  // evolve with RK4 and write to file
  filename_pos = "RK4_pos.txt";
  std::string filename_vel = "RK4_vel.txt";
  ofile_pos.open(filename_pos, std::ofstream::out | std::ofstream::trunc);
  ofile_vel.open(filename_vel, std::ofstream::out | std::ofstream::trunc);

  trap.write_to_file(filename_pos, filename_vel);
  for (int i = 1; i <N; i++)
  {
    trap.evolve_RK4(dt);
    trap.write_to_file(filename_pos, filename_vel);
  }
  ofile_pos.close();
  ofile_vel.close();

  // reset the particles in the trap back to it's initial conditions
  particle_2 = Particle(q,m,r,v);
  trap.update_particle(particle_2,1);
  r = {20,0,20};
  v = {0,25,0};
  particle_1 = Particle(q,m,r,v);
  trap.update_particle(particle_1,0);

  // computing the motion of the particles without Coulomb interactions
  trap.use_Coulomb_interactions(false);
  filename_pos = "RK4_pos_no_inter.txt";
  filename_vel = "RK4_vel_no_inter.txt";
  ofile_pos.open(filename_pos, std::ofstream::out | std::ofstream::trunc);
  ofile_vel.open(filename_vel, std::ofstream::out | std::ofstream::trunc);

  trap.write_to_file(filename_pos, filename_vel);
  for (int i = 1; i <N; i++)
  {
    trap.evolve_RK4(dt);
    trap.write_to_file(filename_pos, filename_vel);
  }
  ofile_pos.close();
  ofile_vel.close();

  // error calculations

  trap = PenningTrap(B0, V0d, d);
  arma::vec N_err = {4000,8000,16000,32000};
  r = {20,0,20};
  v = {0,25,0};
  particle_1 = Particle(q,m,r,v);
  trap.add_particle(particle_1);
  for (int N_i : N_err)
  {
    double h_i = 1.*T/N_i;

    filename_pos = "RK4_error_calc_N_"+std::to_string(N_i)+".txt";
    ofile_pos.open(filename_pos, std::ofstream::out | std::ofstream::trunc);

    trap.write_to_file(filename_pos);
    for (int i = 1; i <N_i; i++)
    {
      trap.evolve_RK4(h_i);
      trap.write_to_file(filename_pos);
    }
    ofile_pos.close();
    particle_1 = Particle(q,m,r,v);
    trap.update_particle(particle_1,0);

    filename_pos = "fEuler_error_calc_N_"+std::to_string(N_i)+".txt";
    ofile_pos.open(filename_pos, std::ofstream::out | std::ofstream::trunc);
    trap.write_to_file(filename_pos);
    for (int i = 1; i <N_i; i++)
    {
      trap.evolve_fEuler(h_i);
      trap.write_to_file(filename_pos);
    }
    ofile_pos.close();
    particle_1 = Particle(q,m,r,v);
    trap.update_particle(particle_1,0);
  }


  // problem 9

  arma::vec f = {0.1,0.4,0.7}; // amplitude
  int N_omega = (2.5-0.2)/domega;
  arma::vec omega_v = arma::linspace(0.2,2.5,N_omega); // angular frequency
  int M = 100; // number of particles

  T = 500;
  dt = 0.1;
  N = T/dt;
  int part_in_trap;

  t = arma::linspace(0, T, N); // makes t values with stepsize dt

  // create a new trap
  PenningTrap trap2 = PenningTrap(B0, V0d, d);
  // set Coulomb interactions and time dependent potential
  trap2.use_Coulomb_interactions(false);
  trap2.time_dependent_potential(true);
  double t_total;
  clock_t t1,t2;

  // create the particles and add them to the trap
  // with random positions and velocities
  for (int i = 0; i < M; i++)
  {
    r = arma::vec(3).randn()*0.1*d;
    v = arma::vec(3).randn()*0.1*d;
    // create a particle of the Particle class
    Particle new_particle = Particle(q, m, r, v);
    // add the created particle into the Penning trap
    trap2.add_particle(new_particle);
  }

  // evolve the particles in the trap for each amplitude and angular frequency
  t1 = clock();
  for (int j = 0; j < f.size(); j++)
  {
    std::string filename = "part_in_trap_f_"+std::to_string(f(j))+".txt";
    std::ofstream ofile;
    ofile.open(filename);
    std::cout << std::endl << std::endl;
    std::cout << "f = " << f(j) << std::endl;

    for (int k = 0; k < omega_v.size(); k++)
    {
      for (int i = 0; i < N; i++)
      {
        trap2.evolve_RK4(dt,t(i),f(j),omega_v(k));
      }
      part_in_trap = trap2.particles_in_trap();
      ofile << std::setw(12) << std::scientific << omega_v(k)
      << std::setw(6) << std::scientific << part_in_trap << std::endl;
      std::cout << "omega_v = " << omega_v(k) << "   " << "Particles remaining :"
      << part_in_trap << std::endl;
      for (int i = 0; i < M; i++)
      {
        r = arma::vec(3).randn()*0.1*d;
        v = arma::vec(3).randn()*0.1*d;
        Particle updated_particle = Particle(q, m, r, v);
        trap2.update_particle(updated_particle, i);
      }
    }        
    ofile.close();
  }
  t2 = clock();
  t_total += t2 - t1;
  double time = ((double) (t_total))/CLOCKS_PER_SEC;
  std::cout << time << std::endl;

  // investigate one of the resonances found from earlier
  trap2.use_Coulomb_interactions(true);
  double f_i = 0.4; // chosen amplitude
  // chosen frequencies
  double omega_start = 1.15;
  double omega_end = 1.55;
  std::string filename = "C_i_resonans_trap.txt";
  std::ofstream ofile;
  ofile.open(filename);
  N_omega = 20;
  omega_v = arma::linspace(omega_start,omega_end,N_omega);
  for (int k = 0; k < N_omega; k++)
  {
    for (int i = 0; i < N; i++)
    {
      trap2.evolve_RK4(dt,t(i),f_i,omega_v(k));
    }
    part_in_trap = trap2.particles_in_trap();
    ofile << std::setw(12) << std::scientific << omega_v(k)
    << std::setw(6) << std::scientific << part_in_trap << std::endl;
    std::cout << "omega_v = " << omega_v(k) << "   " << "Particles remaining :"
    << part_in_trap << std::endl;
    for (int i = 0; i < M; i++)
    {
      r = arma::vec(3).randn()*0.1*d;
      v = arma::vec(3).randn()*0.1*d;
      Particle updated_particle = Particle(q, m, r, v);
      trap2.update_particle(updated_particle, i);
    }       
  }      
  ofile.close();

  // repeat as above, but without Coulomb interactions
  trap2.use_Coulomb_interactions(false);
  filename = "resonans_trap.txt";
  ofile.open(filename);
  omega_v = arma::linspace(omega_start,omega_end,N_omega);
  for (int k = 0; k < N_omega; k++)
  {
    for (int i = 0; i < N; i++)
    {
      trap2.evolve_RK4(dt,t(i),f_i,omega_v(k));
    }
    part_in_trap = trap2.particles_in_trap();
    ofile << std::setw(12) << std::scientific << omega_v(k)
    << std::setw(6) << std::scientific << part_in_trap << std::endl;
    std::cout << "omega_v = " << omega_v(k) << "   " << "Particles remaining :"
    << part_in_trap << std::endl;
    for (int i = 0; i < M; i++)
    {
      r = arma::vec(3).randn()*0.1*d;
      v = arma::vec(3).randn()*0.1*d;
      Particle updated_particle = Particle(q, m, r, v);
      trap2.update_particle(updated_particle, i);
    }       
  }      
  ofile.close();

  return 0;
}
