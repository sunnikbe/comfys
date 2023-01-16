
#include "tridiagonal.hpp"
#include "analytical_eigen.hpp"

// short script that test the create tridiagonal matrix function from
// "tridiagonal.hpp" script and then saves it as a txt file
int main(){
    int N = 6;
    int n = 10;
    double h = std::pow(1./n,2);
    double a = -1;
    double d = 2;

    arma::mat A = create_symmetric_tridiag(N,a,d);    

    // use the eig_sym function for armadillo to find the eigenvalues and
    // eigenevectors for the matrix A
    arma::mat eigvec;
    arma::vec eigval;
    eig_sym(eigval, eigvec, A);

    printf("armadillo:\n");
    printf("Eigenvalues:\n");
    eigval.print(std::cout);
    printf("Eigenvectors\n");
    eigvec.print(std::cout);

    // computes the eigenvectors and eigenvalues analytically
    // for a symmetric tridiagonal matrix and writes them to txt files
    // "eig_val_anal.txt" & "eig_vec_anal.txt"
    compute_analytical_eigen_val_vec(N,d,a);

    return 0;
}