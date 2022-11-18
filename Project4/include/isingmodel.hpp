
#ifndef __isingmodel_hpp__
#define __isingmodel_hpp__

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


class Isingmodel
{
    public:
        arma::mat spin,E_mat,M_mat,EE_mat,MM_mat;
        double T,prob[17];
        int L,E,M,seed;


        std::mt19937 generator;
        std::uniform_int_distribution<int> uniform_int;
        std::uniform_real_distribution<double> uniform_dist;

        // Constructor
        Isingmodel(arma::mat, double, int);

        // Methods
        void print_state();

        void set_spin_state(arma::mat spin_in);

        double compute_prob_factor(int delta_E);

        int compute_total_energy(arma::mat spin);

        int compute_total_magnetization(arma::mat spin);

        int compute_delta_E(int i, int j, arma::mat spin_in);

        void MCMC(int index, int thread_id, int E_thread, int M_thread, int dE,
        arma::mat spin_thread);

        void compute_expected_values(std::string filename, int cycles);

};

#endif