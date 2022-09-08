#include "header.hpp"

// Defining x vector and filling with values between 0 and 1
int N = 100;
arma::vec x = arma::linspace(0.0, 1.0, N + 1);

// declaring x_i for some reason
double x_i;

double u(double x_i) {
  return 1 - (1 - exp(-10)*x_i - exp(-10*x_i));
}

int main(){
  // setting filename for txt file
  std::string filename = "p2.txt";

  // creating and opening output file
  std::ofstream ofile;
  ofile.open(filename);

  // formatting output
  int width = 12, prec = 4; // four decimals and 12 characters

  // function
  double res[N]; //should return our function values

  // Loop over steps
    for (int i = 0; i <= (N); i++)
    {
      res[i] = u(x[i]);
    // Write a line with the current x and y values (nicely formatted) to file
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << res[i]
          << std::endl;
  }

  // closing output file
  ofile.close();
  return 0;
}
