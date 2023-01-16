## Project 3

#### Collaboration by Lars Opgård & Sunniva Kiste Bergan

# Content

The project contains one main c++ file, two c++ header files and one python file
 
    main.cpp
    The main c++ file for this project. Creates particles through the Particle.cpp file and computes
    the motion of the particles in the Penning trap through the PenningTrap.cpp file.
    
    Particle.cpp
    c++ program that contains a Particle class built from the charge, mass, position and velocity of a particle
    
    PenningTrap.cpp
    c++ program that calculates the force that affect a particle inside of a Penning trap and then evolves them using
    Runge-Kutta4 and forward Euler. The class includes functions; write to file, compute with and without
    Coulomb interactions or time dependent electric potential and also count how many particles remains inside the trap.
    
    plots.py
    A python program that plots the results from main.cpp. Also includes functions to compute the analytically solution for 
    the one particle Penning trap system, the relative error between the numerical approximations and analytical solution,
    as well as compute the error convergence rate.
    


# Compile and run code

    main.cpp
    g++ main.cpp -std=c++11 src/*.cpp -I include -o main.exe -larmadillo
    ./main.exe
    
    plots.py
    python3 plots.py
