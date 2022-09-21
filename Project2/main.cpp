
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <assert.h>
#include "tridiagonal.hpp"

double largest_off_diagonal_element(arma::mat A, int& k, int& l);
void jacobi_rotate(arma::mat& A, arma::mat& R);
void jacobi_eigensolver();
int main(){
    // create a symmetric matrix A of size NxN
    int N = 6;
    //arma::mat A = arma::mat(N, N).zeros();
    
    double a = 1;
    double d = 2;
    arma::mat A = create_symmetric_tridiag(N,a,d);  

    // test the largest off diagonal element function
    arma::mat B = arma::mat(4,4).zeros();
    B.diag().ones();
    B(1,2) = -0.7;
    B(2,1) = B(1,2);
    B(0,3) = 0.5;
    B(3,0) = B(0,3);

    double max_val_test;
    int k,l;
    max_val_test = largest_off_diagonal_element(B,k,l);

    printf("max value = %4.2f \n", max_val_test);
    printf("k = %i \n", k);
    printf("l = %i \n", l);

    // create matrix R^(1) = I, R^(m) = S_m
    arma::mat R = arma::mat(N,N,arma::fill::eye);
    jacobi_rotate(A, R);
    A.save("diag_mat.txt", arma::raw_ascii);

    std::string filename = "eigenvectors.txt";
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::setw(-1) << "#" << std::setw(7) << "v_1" << std::endl;
    ofile.close();
    R.save(filename, arma::raw_ascii);

    arma::vec eig_val = A.diag();
    filename = "eigenvalues.txt";
    ofile.open(filename);
    ofile << std::setw(-1) << "#" << std::setw(13) << "eigenvalues" << std::endl;
    for (int i = 0; i < N; i++){
        ofile << std::setw(15) << std::scientific << eig_val(i) << std::endl;
    }

    return 0;
}

double largest_off_diagonal_element(arma::mat A, int& k, int& l){
    int n = A.n_rows;

    // test that A is a square matrix, gives error if not
    assert (A.is_square() != 0);

    double max_val;
    k = 0;
    l = 1;
    for (int j = 1; j < n; j++){
        for (int i = 0; i < j; i++){
            if (std::abs(A(i,j)) > max_val){
                max_val = std::abs(A(i,j));
                // set the new indeces to k and l
                k = i;
                l = j;
            }
        }
    }

    return max_val;
}

// solve the eigen equation using Jacobi rotation method
void jacobi_rotate(arma::mat& A, arma::mat& R){
    double epsilon = 1e-8; // toleranse = epsilon
    // declare variables used in the algorithm
    double tau, tan, cos, sin, max_val;
    int N = A.n_rows;
    int k,l;
    // start with finding the largest off-diagonal value in A
    max_val = largest_off_diagonal_element(A,k,l);

    int num_iter = 0; // count the number of iterations done
    while(max_val > epsilon){
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0){
            tan = 1/(tau + std::sqrt(1 + std::pow(tau,2)));
        }
        else{
            tan = 1/(-tau + std::sqrt(1 + std::pow(tau,2)));
        }
        cos = 1/(1 + std::pow(tan,2));
        sin = cos*tan;

        // transforming the current A matrix
        // by first updating the elements with indeces k and l
        A(k,k) = A(k,k)*std::pow(cos,2) - 2*A(k,l)*cos*sin + A(l,l)*std::pow(sin,2);
        A(l,l) = A(l,l)*std::pow(cos,2) + 2*A(k,l)*cos*sin + A(k,k)*std::pow(sin,2);
        A(k,l) = 0;
        A(l,k) = 0;

        // then updating the rest of the elements where i != k,l along the k,l rows and columns
        double a_ik, a_il;
        for (int i = 0; (i != k, i != l, i < N); i++){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik*cos - a_il*sin;
            A(k,i) = A(i,k);
            A(i,l) = a_il*cos + a_ik*sin;
            A(l,i) = A(i,l);
        }

        // update the overall rotation matrix
        double r_ik;
        for (int i = 0; i < N; i++){
            r_ik = R(i,k);
            R(i,k) = R(i,k)*cos - R(l,l)*sin;
            R(i,l) = R(i,l)*cos + R(i,k)*sin;
        }

        // find the next off-diagonal element and continue until the tolerance is reached
        max_val = largest_off_diagonal_element(A,k,l);
        num_iter += 1;
    }
    
    printf("number of iterations = %i \n", num_iter);
}