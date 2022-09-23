
#ifndef __tridiagonal_hpp__
#define __tridiagonal_hpp__

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <assert.h>

// function creates the tridiagonal matrix A(a,d,e) of size NxN 
// from three vectors using the diag() function from armadillo
arma::mat create_tridiag(int N, arma::vec a_vec, arma::vec d_vec, arma::vec e_vec);

// takes floats and create vectors for the tridiagonal
arma::mat create_tridiag(int N, double a, double d, double e);

// creates a symmetric tridiagonal matrix A(a,d,a) of size NxN
arma::mat create_symmetric_tridiag(int N, double a, double d);

#endif