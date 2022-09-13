#include "header.hpp"

//Special algorithm (a lot is from my general algorithm code)
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
arma::vec gt = arma::vec(m);


int main(){
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

  // Writing values to txtfile
  // This part is copy pasted from p2.cpp and changed

  // setting filename for txt file
  std::string filename = "p9.txt";

  // creating and opening output file
  std::ofstream ofile;
  ofile.open(filename);

  // formatting output
  int width = 12, prec = 4; // four decimals and 12 characters

  // Loop over steps
    for (int i = 0; i < m; i++){
      // Write a line with the current x and y values (nicely formatted) to file
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
            << std::setw(width) << std::setprecision(prec) << std::scientific << x[i]
            << std::endl;
  }

  // closing output file
  ofile.close();
  return 0;
}
