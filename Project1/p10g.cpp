#include "header.hpp"
#include <time.h>

// Timing general alg. p7.cpp

//General algorithm
double n_steps = 10.; // Is to be changed

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
//tilde vectors
arma::vec bt = arma::vec(m);
arma::vec gt = arma::vec(m);

int main ()
{
  // Start measuring time
  clock_t t1 = clock();





  //Initial values
  bt[0] = b;
  gt[0] = hh*f(x[0]);

  // Forward sub
  for (int i = 1; i <= (m - 1); i++){
    bt[i] = b - (a/(bt[i - 1]))*c;
    gt[i] = hh*f(x[i]) - (a/(bt[i - 1]))*gt[i - 1];
  }


  // Back sub
  for (int i = (m - 2); i >= 1; i--){
    v[i] = (gt[i] - (c*(v[i + 1])))/(bt[i]);
  }

  //Last element of v - vec
  v[m] = (gt[m])/(bt[m]);





  // Stop measuring time
  clock_t t2 = clock();

  // Calculate the elapsed time.
  double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

  // Write to file for calc. later
  // setting filename for txt file
  std::string filename = "p10.txt";

  // creating and opening output file
  std::ofstream ofile;
  ofile.open(filename); //, std::ofstream::app);

  // formatting output
  int width = 12, prec = 4; // four decimals and 12 characters

  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << n_steps
        << std::setw(width) << std::setprecision(prec) << std::scientific << duration_seconds
        << std::endl;
  }

  // closing output file
  ofile.close();
  return 0;
}
