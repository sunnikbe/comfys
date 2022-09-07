// minimal example file

#include <iostream>
#include <string>
#include <fstream>

int main()
{
  // Set a filename
  std::string filename = "output.txt";

  // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect
  // it to our filename.
  std::ofstream ofile;
  ofile.open(filename);

  // Write: ofile.open(filename, std::ofstream::app); if you want to append to
  // a file and not overwrite

  // send some text top this output file
  ofile << "Some output text" << std::endl;

  // close the oputput file
  ofile.close();

  // All is well. Exit program with return code 0.
  return 0;
}


// For the terminal/README file:

//g++ outputs.cpp -o outputs.exe
//./outputs.exe
