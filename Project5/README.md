## Project 5

#### Collaboration by Lars Opgård & Sunniva Kiste Bergan

# Content

The project contains three c++ files, one python file and one mp4 file

  main.cpp \
  The main c++ program that uses the simulate_slit_experiment.cpp and simulates the \
  wave packet thorugh a slit experiment for none or more slit openings
  
  simulate_slit_experiment.cpp \
  c++ program that includes the class that simulates the slit experiment for a given \
  initial state and slit experiment
  
  simulate_slit_experiment.hpp \
  c++ header file for the simulate_slit_experiment.cpp program
  
  plots.py \
  python program that plots the results given by the main.cpp program. Uses terminal input \
  for number if slit openings, only used for naming files
  
  
# Compile and run

  main.cpp \
  g++ main.cpp -std=c++11 src/*.cpp -I include -o main.exe -larmadillo \
  ./main.exe "number_of_slits" "potential strenght" "sigma_y" "T"
  
  where the sigma_y is the width of the wave packet in the y direction, \
  T is the end time for the simulation

  plots.py \
  python3 plots.py
  
