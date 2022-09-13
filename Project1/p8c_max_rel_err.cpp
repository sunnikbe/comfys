#include "header.hpp"

//General algorithm
double n_steps = pow(10., 7.); // Is to be changed

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

// True value eq.
double u(double x_j) {
  return 1 - (1 - exp(-10))*x_j - exp(-10*x_j);
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


int main(){
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

  arma::vec real_err = arma::vec(m);
  for (int j = 1; j <= (m - 1); j++){
    real_err[j] = abs((v[j] - u(x[j]))/(u(x[j])));
  }

  double max_rel_err = max(real_err);

  // Writing values to txtfile
  // This part is copy pasted from p2.cpp and changed

  // setting filename for txt file
  std::string filename = "n_steps_max_rel_err.txt";

  // creating and opening output file
  std::ofstream ofile; // add under when ready
  ofile.open(filename, std::ofstream::app);

  // formatting output
  int width = 20, prec = 10; // four decimals and 12 characters

  // Write n_steps and max_rel_err to file
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << n_steps
        << std::setw(width) << std::setprecision(prec) << std::scientific << max_rel_err
        << std::endl;

  // closing output file
  ofile.close();
  return 0;
}
