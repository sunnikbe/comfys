// creating a function

// type function_name(type1 input1, ..., typeN inputN){ //declaration
    //Code block or body (i.e definition)
//}

// f(x) = x*x function code

#include <iostream>

double f(double x); // decleration

int main() {
  double x = 2;
  double res = f(x); // should return 2*2 = 4.
  std::cout << "f(x) = " << res << "\n"; // prints f(x) = 4
  return 0;
}

double f(double x){
  // Definition of function. It squares stuff.
  return x*x;
}
