
#include "analytical_eigen.hpp"

// computes the eigenvalues and eigenvectors only for a symmetric tridiagonal matrix
void compute_analytical_eigen_val_vec(int N, double d, double a){
    arma::vec eig_val = arma::vec(N).zeros();
    arma::mat eig_vec = arma::mat(N,N).zeros();
    for (int j = 1; j <= N; j++){
        eig_val(j-1) = d + 2*a*std::cos(j*M_PI/(N+1));
        for (int i = 1; i <= N; i++){
            eig_vec(i-1,j-1) = std::sin(i*j*M_PI/(N+1));
        }
    }
    eig_vec = arma::normalise(eig_vec);

    printf("Analytical\n");
    printf("Eigenvalues:\n");
    printf("Eigenvectors:\n");
    eig_val.print(std::cout);
    eig_vec.print(std::cout);
}