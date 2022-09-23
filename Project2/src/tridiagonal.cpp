
#include "tridiagonal.hpp"

// function creates the tridiagonal matrix A(a,d,e) of size NxN 
// from three vectors using the diag() function from armadillo
arma::mat create_tridiag(int N, arma::vec a_vec, arma::vec d_vec, arma::vec e_vec){
    arma::mat A = arma::mat(N,N).zeros();

    A.diag() = d_vec; // main diagonal
    A.diag(-1) = a_vec; // sub diagonal
    A.diag(1) = e_vec; // super diagonal

    return A;
}

// takes floats and create vectors for the tridiagonal
arma::mat create_tridiag(int N, double a, double d, double e){
    arma::vec a_vec = arma::vec(N-1).ones() * a;
    arma::vec d_vec = arma::vec(N).ones() * d;
    arma::vec e_vec = arma::vec(N-1).ones() * e;

    // calls the function above that creates the matrix
    return create_tridiag(N, a_vec, d_vec, e_vec);
}

// creates a symmetric tridiagonal matrix A(a,d,a) of size NxN
arma::mat create_symmetric_tridiag(int N, double a, double d){
    double e = a;

    return create_tridiag(N, a, d, e);
}
