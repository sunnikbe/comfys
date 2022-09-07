// example code
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{

  //
  // Check number of command-line arguments
  //

  if (argc != 4)  // Expect 3 command-line arguments
  {
    // Get the name of the executable file
    std::string executable_name = argv[0];

    std::cerr << "Error: Wrong number of input arguments." << std::endl;
    std::cerr << "Usage: " << executable_name << " <some integer> <some double> <some string>" << std::endl;

    // Exit program with non-zero return code to indicate a problem
    return 1;
  }

  //
  // Save the command-line arguments as different variables
  //

  int some_integer = atoi(argv[1]);    // atoi converts argv[2] to an integer
  double some_double = atof(argv[2]);  // atof converts argv[3] to a floating-point number
  std::string some_string = argv[3];   // no type conversion needed


  //
  // Some dummy code that uses the input arguments
  //

  std::cout << "Got the following input: "
            << "some_integer: " << some_integer << ", "
            << "some_double: " << some_double << ", "
            << "some_string: " << some_string
            << std::endl;


  // All is well. Exit program with return code 0.
  return 0;
}

// Terminal:
// Sunnivas-MacBook-Pro:comfys sunnivakbergan$ ./inputs.exe 100 3.14 hello

// Produced the output:
// Got the following input: some_integer: 100, some_double: 3.14, some_string: hello
