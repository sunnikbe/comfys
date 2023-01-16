## Project 2

#### Collaboration by Lars Opgård & Sunniva Kiste Bergan

### Content

The project contains two header files, two c++ files and two python files 

The two header files are  

    - tridiagonal.hpp  
    c++ code that creates a NxN symmetric tridiagonal matrix from the elements d
    on the diagonal and a on the sub- superdiagonal  
    
    - analytical_eigen.hpp  
    c++ code that computes the eigenvalues and eigenvectors analytically
    for a symmetric tridiagonal matrix and then print then out to the terminal  
    
The two c++ files are  

    - tridiag.cpp  
    c++ code that uses the tridiagonal.hpp to create a simple symmetric tridiagonal matrix
    and compares the eigenvalues and eigenvectors given by the armadillo eig_sym function
    compared to the analytically computed eigenvalues and eigenvectors from analytical_eigen.hpp  
    
    - main.cpp  
    The main c++ program that solves the eigenfunction using the Jacobi rotate algorithm on a NxN symmetric tridiagonal matrix  
    
The python files are  

    - Nvsnum_iter.py  
    Python code that plots the number of iterations needed to compute the
    Jacobi rotate method as a function of the matrix size N, also creates a polyfit
    for the data  
    
    -buckling_beam.py  
    Plots the eigenvectors from the Jacobi rotate method for
    the three smallest eigenvalues as a function of position x  


### Compile and run code

    main.cpp
    g++ main.cpp -std=c++11 src/*.cpp -I include -o main.exe -larmadillo
    ./main.exe

    tridiag.cpp
    g++ tridiag.cpp -std=c++11 src/*.cpp -I include -o tridiag.exe -larmadillo
    ./tridiag.exe
