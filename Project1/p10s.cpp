#include "header.hpp"
#include <time.h>

// Timing special alg. p9.cpp

//General algorithm
double n_steps = pow(10.,6.); // Is to be changed

double m = n_steps + 1; //m from problem 5
double h = 1/n_steps; // stepsize
double hh = h*h;
double n = m - 2; //matrix A is nxn

// x - vector
arma::vec x = arma::linspace(0.0, 1.0, m);

// Poisson f(x):
double f(double x_i) {
  return 100.0*std::exp(-10.0*x_i);
}

//Defining a, b, c in matrix A
double a = -1.0;
double b = 2.0;
double c = -1.0;

//Defining v and g vectors
arma::vec v = arma::vec(m);
arma::vec g = arma::vec(m);
arma::vec gt = arma::vec(m);

int main ()
{
  // Start measuring time
  clock_t t1 = clock();





  //Initial conditions
  gt[0] = hh*f(x[0]);
  v[m] = 0.0;

  // Forward sub
  double d = a/b;
  for (int i = 1; i <= (m - 1); i++){
    gt[i] = hh*f(x[i]) - d*gt[i - 1];
  }

  // Back sub
  for (int i = (m - 2); i >= 1; i--){
    v[i] = (gt[i] - c*(v[i + 1]))/b;
  }





  // Stop measuring time
  clock_t t2 = clock();

  // Calculate the elapsed time.
  double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

  // Write to file for calc. later
  // setting filename for txt file
  std::string filename = "p10s.txt";

  // creating and opening output file
  std::ofstream ofile;
  ofile.open(filename, std::ofstream::app);

  // formatting output
  int width = 20, prec = 10; // 10 decimals and 20 characters

  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << n_steps
        << std::setw(width) << std::setprecision(prec) << std::scientific << duration_seconds
        << std::endl;

  // closing output file
  ofile.close();
  return 0;
}
