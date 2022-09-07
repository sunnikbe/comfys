// you need to include the header file:
#include <armadillo>

// Just some bs to test
#include <iostream>

int main() {
  std::cout << "It works\n";
  return 0;
}

// When compiling you have to add the compiler flag:
// -std=c++11
// during compilation:
// g++ -c intro_armadillo.cpp -std=c++11

// When linking, compiler flag:
// -larmadillo
// is added during linking:
// g++ intro_armadillo.o -o intro_armadillo.exe -larmadillo

// put together:
// g++ intro_armadillo.cpp -std=c++11 -o intro_armadillo.exe -larmadillo
