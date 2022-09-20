
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <assert.h>
#include "tridiagonal.hpp"
#include "analytical_eigen.hpp"

// short script that test the create tridiagonal matrix function from
// "tridiagonal.hpp" script and then saves it as a txt file
int main(){
    int N = 6;
    double a = 1;
    double d = 2;

    arma::mat A = create_symmetric_tridiag(N,a,d);    

    A.save("Output.txt", arma::raw_ascii);

    // use the eig_sym function for armadillo to find the eigenvalues and
    // eigenevectors for the matrix A
    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, A);

    eigvec.save("eigvec.txt", arma::raw_ascii);
    eigval.save("eigval.txt", arma::raw_ascii);

    compute_analytical_eigen_val_vec(N,d,a);

    return 0;
}