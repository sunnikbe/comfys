
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
        bool ordered;

        std::mt19937 generator;
        std::uniform_int_distribution<int> uniform_int,spin_dist;
        std::uniform_real_distribution<double> uniform_dist;

        // Constructor
        Isingmodel(double, int, int, bool ordered = true);

        // Methods
        int compute_total_energy(arma::mat spin);

        int compute_total_magnetization(arma::mat spin);

        void mcmc(std::string filename, int cycles, int burn_in = 0, bool mat_form = true);
};

#endif
