
#include "simulate_slit_experiment.hpp"

using namespace std;
using namespace arma;

Slit_simulation::Slit_simulation(int M_in, double h_in, double dt_in, double T_in,
                                 double v0_in, int number_of_slits_in)
{
    M = M_in; // number of spatial steps
    N = (M-2)*(M-2); // size of the u vector, and A and B have size NxN
    // wihtout the boundary conditions

    h = h_in; // spatial stepsize
    dt = dt_in; // time stepszie
    T = T_in; // endtime of the simulation
    n = T/dt; // number of timesteps

    // spatial vectors wiht boundary conditions
    x = linspace(0,1,M);
    y = linspace(0,1,M);
    u.zeros(N); // u vector where u_k = u_ij, find k using the find_k(int i, intj) function
    // a and b vector used in creating the A and B matrices
    a.zeros(N);
    b.zeros(N);

    // values used for creating the slits
    thickness = 0.02; // wall thickness from the center of the wall in x direction
    centre = 0.5; // position of the center in the x direction
    length = 0.05; // length of the wall seperating slits, only used for 2 or more slits
    slit_size = 0.05; // lenght of the opening
    // the slits will be symmetric around y = 0.5
    v0 = v0_in; // strength of the potential wall
    number_of_slits = number_of_slits_in;
    V = mat(M,M,fill::zeros); // matrix mapping the potential wall
    // if the potential is set to 0, the program will run the compute_potential() function
    // if however v0 != 0 and number of slits is set to 0, the program will make a holeless wall
    // if you would ever need that for some reason
    if (v0 != 0)
    {
        compute_potential(); // map the potential wall
    }
}

sp_cx_mat Slit_simulation::create_diagonal_matrix(cx_vec d, complex<double> r)
{
    // create a vector that becomes the main diagonal of A and B matrix
    cx_vec r_vec(N-1, fill::value(r));
    for (int i = M-3; i < N-1; i += M-2)
    {
        r_vec(i) = 0;
    }
    // create the matrix
    sp_cx_mat diag_mat(N,N);
    diag_mat.diag(0) = d;
    diag_mat.diag(M-2) += r;
    diag_mat.diag(-(M-2)) += r;
    diag_mat.diag(1) = r_vec;
    diag_mat.diag(-1) = r_vec;
    return diag_mat;
}

int Slit_simulation::find_k(int i, int j)
{
    // return index k for a given index i and j, basically returns the index k for an
    // element c_k from the flattened matrix with position c_ij
    return i + j*(M-2);
}

void Slit_simulation::compute_potential()
{
    // set up the wall at given x position centre
    double b0 = centre - thickness; // start of wall
    double b1 = centre + thickness; // end of wall

    // find the indices for the start of the lowest slit to the end of the highest slit
    // then compute the index difference between the length of the slit openings
    int n_slit0 = (0.5 - (2*number_of_slits - 1)/2*slit_size)/1.0*M;
    int n_slit1 = (0.5 + (2*number_of_slits - 1)/2*slit_size)/1.0*M;
    int n_dslit = slit_size/h;
    vec v_i(M); // declare a column vector equal to the size of a column of V
    // for each index where the wall is present in the x direction
    for (int i = b0/1.0*M; i < b1/1.0*M; i++)
    {
        v_i.fill(v0); // fill a column vector with the potential strength
        // for each slit used in the simulation
        for (int j = 0; j < number_of_slits; j++)
        {
            int n0 = n_slit0 + 2*j*n_dslit; // start of the slit
            int n1 = n_slit0 + n_dslit + 2*j*n_dslit - 1; // end of the slit
            v_i(span(n0,n1)).zeros(); // set the span between the slit opening to zero
        }
        V.col(i) = v_i; // set the column of V equal to the potential vector
    }
    // save the potential map for later if necessary
    V.save("Potential_map.txt",raw_ascii);
}

void Slit_simulation::find_A_and_B()
{
    // compute r = i*dt/2h^2
    r = complex<double>(0,1.0)*dt/2.0/h/h;
    // compute a and b values for A and B respectably
    for (int j = 0; j < M-2; j++)
    {
        for (int i = 0; i < M-2; i++)
        {
            int k = find_k(i,j);
            a(k) = 1.0 + 4.0*r + complex<double>(0,1.0)*dt/2.0*V(j+1,i+1);
            b(k) = 1.0 - 4.0*r - complex<double>(0,1.0)*dt/2.0*V(j+1,i+1);
        }
    }
    // create A and B matrix
    A = create_diagonal_matrix(a,-r);
    B = create_diagonal_matrix(b,r);
}

// function that updates the U matrix for the current time step
void Slit_simulation::update_U(int slice)
{
    for (int j = 1; j < M-1; j++)
    {
        for (int i = 1; i < M-1; i++)
        {
            int k = find_k(i-1,j-1);
            U(j,i,slice) = u(k);
        }
    }
}

// compute the intital state using a Gaussian wave packet
void Slit_simulation::initial_state(double xc, double yc, double px, double py, double sigmax, double sigmay)
{
    for (int j = 0; j < M-2; j++)
    {
        for (int i = 0; i < M-2; i++)
        {
            int k = find_k(i,j);
            u(k) = exp(-((x(i)-xc)*(x(i)-xc))/(2*sigmax*sigmax) - ((y(j)-yc)*(y(j)-yc))/(2*sigmay*sigmay)
            + complex<double>(0,1.0)*px*((x(i)-xc)) + complex<double>(0,1.0)*py*(y(j)-yc));
        }
    }
    // normalize the initial state
    complex<double> norm = cdot(u,u);
    u = u/sqrt(norm);
    U.zeros(M,M,n+1); // set all the elements in U to zero
    update_U(0); // update U with the initial state
}

void Slit_simulation::evolve_next_time_step()
{
    // compute A and B matrices
    find_A_and_B();
    // evolve the state for each time step in two steps
    for (int i = 0; i <= n; i++)
    {
        // first compute the b vector using matrix vector multiplication between B and u^n
        cx_vec b_vec = B*u;
        // solve the equation Au^(n+1) = b for the next time step of the state u^(n+1)
        spsolve(u,A,b_vec);
        // update U
        update_U(i);
    }
    // save the states as a binary file
    U.save("Output.bin",arma_binary);
}
