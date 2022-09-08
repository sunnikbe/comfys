#include "header.hpp"

// Defining x vector and filling with values between 0 and 1
int N = 100;
arma::vec x = arma::linspace(0.0, 1.0, N);

// declaring x_i for some reason
double x_i;

double u(double x_i) {
  return 1 - (1 - exp(-10)*x_i - exp(-10*x_i));
}

int main(){
  // setting filename for txt file
  std::string filename = "p1p2data.txt";

  // creating and opening output file
  std::ofstream ofile;
  ofile.open(filename);

  // formatting output
  int width = 12, prec = 4; // four decimals and 12 characters

  // function
  double res = u(x_i); //should return our function values

  // Loop over steps
    for (int i = 0; i <= (N - 1); i++)
    {
    // Write a line with the current x and y values (nicely formatted) to file
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << double (x_i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << res
          << std::endl;

    // Function
    double x_i = x(i);
  }

  // closing output file
  ofile.close();
  return 0;
}
