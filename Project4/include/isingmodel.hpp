
#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>


class Isingmodel
{
  public:
    arma::mat spin;
    double T,dE0,dE4p,dE4m,dE8p,dE8m;
    int L,E_before,M;
    std::mt19937 generator;
    std::uniform_int_distribution<int> uniform_int;
    std::uniform_real_distribution<double> uniform_dist;

    // Constructor
    Isingmodel(arma::mat, double, int);

    // Methods
    void print_state();

    double compute_prob_factor(int delta_E);

    int compute_total_energy(arma::mat spin);

    int compute_total_magnetization(arma::mat spin);

    void spin_candidate();

    void MarkovChainMonteCarlo(std::string filename, int iterations);

};

#endif