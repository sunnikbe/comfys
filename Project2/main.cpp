
#include "tridiagonal.hpp"

double largest_off_diagonal_element(arma::mat A, int& k, int& l);
void jacobi_rotate(arma::mat& A, arma::mat& R,int& num_iter);
void jacobi_eigensolver(arma::mat &A, arma::mat& R, arma::mat& eigenvectors, arma::vec& eigenvalues, int& max_iter, int& num_iter, double tol);
int main(){
    std::string filename;
    std::ofstream ofile;
    arma::mat eigvecs;
    arma::vec eigvals;
    int max_iter = 100000;
    double epsilon = 1e-8;
    int N, n, num_iter, k, l;
    double h, a, d;

    // Discretization with n = 10 steps and n = 100 steps

    n = 100; // number of steps
    N = n - 1; // interior points
    h = 1./n; // stepsize

    // create a symmetric matrix A of size NxN
    //arma::mat A = arma::mat(N, N).zeros();

    double h_sq = std::pow(h, 2);

    a = -1/h_sq;
    d = 2/h_sq;
    arma::mat A = create_symmetric_tridiag(N,a,d);
    arma::mat R = arma::mat(N, N, arma::fill::eye);
    jacobi_eigensolver(A, R, eigvecs, eigvals, max_iter, num_iter, epsilon);
    eigvecs.save("eigenvectors.txt", arma::raw_ascii);
    eigvals.save("eigenvalues.txt", arma::raw_ascii);

    // // Problem 5a,b) running jacobi_rotate with diagonal matrix and dense
    // // Writing to file
    // filename = "dense.txt"; // For 5a) "Nvsnum_iter.txt";
    // ofile.open(filename);
    //
    // arma::vec N_vals = {25, 50, 75, 100, 125, 150, 175, 200};
    // for (int i: N_vals){
    //     // problem 5a)
    //     // arma::mat A = create_symmetric_tridiag(i, a, d);
    //
    //     // problem 5b)
    //     arma::mat A = arma::mat(i, i).randn();
    //     // Symmetrize the matrix by reflecting the upper triangle to lower triangle
    //     A = arma::symmatu(A); // only for 5b
    //     arma::mat R = arma::mat(i, i, arma::fill::eye);
    //     jacobi_eigensolver(A,R,eigvecs,eigvals,max_iter,num_iter,epsilon);
    //
    //     ofile << std::setw(3) << i
    //     << std::setw(20) << std::scientific << num_iter
    //     << std::endl;
    // }
    // ofile.close();


    // Problem 3b) test the largest off diagonal element function
    arma::mat B = arma::mat(4,4).zeros();
    B.diag().ones();
    B(1,2) = -0.7;
    B(2,1) = B(1,2);
    B(0,3) = 0.5;
    B(3,0) = B(0,3);

    double max_val_test;
    max_val_test = largest_off_diagonal_element(B,k,l);

    // print the maximum value and indices k and l
    printf("max value = %4.2f \n", max_val_test);
    printf("k = %i \n", k);
    printf("l = %i \n", l);

    // // problem 4
    // // create matrix R^(1) = I, R^(m) = S_m
    // N = 6;
    // arma::mat R = arma::mat(N,N,arma::fill::eye);
    // A = create_symmetric_tridiag(N,-1,2);
    // jacobi_eigensolver(A,R,eigvecs,eigvals,max_iter,num_iter,epsilon);
    //
    // printf("Eigenvalues:\n");
    // eigvals.print(std::cout);
    // printf("Eigenvectors as columns:\n");
    // eigvecs.print(std::cout);

    return 0;
}

double largest_off_diagonal_element(arma::mat A, int& k, int& l){
    int n = A.n_rows;

    // test that A is a square matrix, gives error if not
    assert (A.is_square() != 0);

    double max_val = 0;
    k = 0;
    l = 1;
    for (int j = 1; j < n; j++){
        for (int i = 0; i < j; i++){
            if (std::abs(A(i,j)) > max_val){
                max_val = std::abs(A(i,j));
                // set the new indices to k and l
                k = i;
                l = j;
            }
        }
    }

    return max_val;
}


// solve the eigen equation using Jacobi rotation method
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    // declare variables used in the algorithm
    double tau, tan, cos, sin;
    int N = A.n_rows;

    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if (tau > 0){
        tan = 1./(tau + std::sqrt(1 + std::pow(tau,2)));
    }
    else{
        tan = -1./(-tau + std::sqrt(1 + std::pow(tau,2)));
    }
    cos = 1./std::sqrt(1 + std::pow(tan,2));
    sin = cos*tan;

    // transforming the current A matrix
    // by first updating the elements with indeces k and l
    double a_kk = A(k,k);
    double a_ll = A(l,l);
    A(k,k) = a_kk*std::pow(cos,2) - 2*A(k,l)*cos*sin + a_ll*std::pow(sin,2);
    A(l,l) = a_ll*std::pow(cos,2) + 2*A(k,l)*cos*sin + a_kk*std::pow(sin,2);
    A(k,l) = 0;
    A(l,k) = 0;

    // then updating the rest of the elements where i != k,l along the k,l rows and columns
    double a_ik, a_il, r_ik, r_il;
    for (int i = 0; i < N; i++){
        if (i != k && i != l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik*cos - a_il*sin;
            A(k,i) = A(i,k);
            A(i,l) = a_il*cos + a_ik*sin;
            A(l,i) = A(i,l);
        }

        // update the overall rotation matrix
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = r_ik*cos - r_il*sin;
        R(i,l) = r_il*cos + r_ik*sin;
    }
}


void jacobi_eigensolver(arma::mat &A, arma::mat& R, arma::mat& eigenvectors, arma::vec& eigenvalues, int& max_iter, int& num_iter, double tol){

    // start with finding the largest off-diagonal value in A
    int k = 0;
    int l = 1;
    double max_val;
    max_val = largest_off_diagonal_element(A,k,l);

    num_iter = 0; // count the number of iterations done
    while(max_val > tol && num_iter < max_iter){
        jacobi_rotate(A, R, k, l);
        num_iter++;
        // find the next off-diagonal element and continue until the tolerance is reached
        max_val = largest_off_diagonal_element(A, k, l);
    }

    printf("number of iterations = %i \n", num_iter);

    eigenvalues = A.diag();
    arma::uvec sorted_indices = arma::sort_index(eigenvalues);
    eigenvalues = sort(eigenvalues);
    eigenvectors = R;
    int N = R.n_rows;
    for (int i = 0; i < N; i++){
        eigenvectors.col(i) = R.col(sorted_indices(i));
    }
}
