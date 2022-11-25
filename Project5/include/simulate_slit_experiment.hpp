
#ifndef __simulate_slit_experiment_hpp__
#define __simulate_slit_experiment_hpp__

#include "omp.h"
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
#include <chrono>
#include <complex>

class Slit_simulation
{
    public:
        int M,N;
        double h,dt;
        arma::cx_mat A,B;
        arma::mat V;
        arma::cx_vec u;
        arma::vec x,y;

        //constructor
        Slit_simulation(int M_in, double h_in, double dt_in);

        //methods
        arma::cx_mat create_diagonal_matrix(arma::cx_vec d, std::complex<double> r);

        int find_k(int i, int j);

        void compute_potential(double V0, double thickness, double centre, double lenght, double slit_size, int number_of_slits);

        void find_A_and_B();

        void initial_state(double xc, double yc, double px, double py, double sigmax, double sigmay);

        void evolve_next_time_step();
};

#endif