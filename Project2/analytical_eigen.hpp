
#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <assert.h>

void compute_analytical_eigen_val_vec(int N, double d, double a){
    arma::vec eig_val = arma::vec(N).zeros();
    arma::mat eig_vec = arma::mat(N,N).zeros();
    for (int j = 1; j <= N; j++){
        eig_val(j-1) = d + 2*a*std::cos(j*M_PI/(N+1));
        for (int i = 1; i <= N; i++){
            eig_vec(i-1,j-1) = std::sin(i*j*M_PI/(N+1));
        }
    }

    eig_val.save("eig_val_anal.txt", arma::raw_ascii);
    eig_vec.save("eig_vec_anal.txt", arma::raw_ascii);
}