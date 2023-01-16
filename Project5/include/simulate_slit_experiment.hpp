
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
        int M,N,n,number_of_slits;
        double h,dt,T,v0,thickness,centre,length,slit_size;
        std::complex<double> r;
        arma::sp_cx_mat A,B;
        arma::cx_vec a;
        arma::cx_vec b;
        arma::mat V;
        arma::cx_vec u;
        arma::vec x,y;
        arma::cx_cube U;

        //constructor
        Slit_simulation(int M_in, double h_in, double dt_in, double T_in,
                        double v0_in, int number_of_slits_in);

        //methods
        arma::sp_cx_mat create_diagonal_matrix(arma::cx_vec d, std::complex<double> r);

        int find_k(int i, int j);

        void compute_potential();

        void find_A_and_B();

        void update_U(int slice);

        void initial_state(double xc, double yc, double px, double py, double sigmax, double sigmay);

        void evolve_next_time_step();
};

#endif