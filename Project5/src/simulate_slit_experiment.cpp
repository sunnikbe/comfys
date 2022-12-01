
#include "simulate_slit_experiment.hpp"

using namespace std;
using namespace arma;

Slit_simulation::Slit_simulation(int M_in, double h_in, double dt_in, double T_in)
{
    M = M_in;
    N = (M-2)*(M-2);
    h = h_in;
    dt = dt_in;
    T = T_in;
    n = T/dt;
    x = linspace(0,1,M);
    y = linspace(0,1,M);
    u.zeros(N);
    a.zeros(N);
    b.zeros(N);
}

void Slit_simulation::print_u()
{
    u.print(cout);
    cout << endl;
}

void Slit_simulation::print_U()
{
    U.print(cout);
    cout << endl;
}

void Slit_simulation::print_norm()
{
    cout << sum(cdot(conj(u),u)) << endl;
    cout << endl;
}

cx_mat Slit_simulation::create_diagonal_matrix(cx_vec d, complex<double> r)
{
    cx_vec r_vec(N-1, fill::value(r));
    for (int i = 2; i < N-1; i += 3)
    {
        r_vec(i) = 0;
    }
    cx_mat diag_mat = diagmat(d);
    diag_mat.diag(3) += r;
    diag_mat.diag(-3) += r;
    diag_mat.diag(1) = r_vec;
    diag_mat.diag(-1) = r_vec;
    return diag_mat;
}

int Slit_simulation::find_k(int i, int j)
{
    return i + j*(M-2);
}

void Slit_simulation::compute_potential(double V0, double thickness, double centre, double length, double slit_size, int number_of_slits)
{
    double b0 = centre - thickness;
    double b1 = centre + thickness;
    V = mat(M,M, fill::zeros);
    if (number_of_slits == 1)
    {
        double s0 = 0.5 + slit_size/2;
        double s1 = 0.5 - slit_size/2;
        for (int i = 0; i < M; i++)
        {
            if (x(i) >= b0 && x(i) <= b1)
            {
                for (int j = 0; j < M; j++)
                {
                    if (y(j) >= s0 || y(j) <= s1)
                    {
                        V(j,i) = V0;
                    }
                }
            }
        }
    }
    if (number_of_slits == 2)
    {
        double s0 = 0.5 + slit_size/2;
        double s1 = 0.5 - slit_size/2;
        double s2 = 0.5 + 3.0*slit_size/2;
        double s3 = 0.5 - 3.0*slit_size/2;
        for (int i = 0; i < M; i++)
        {
            if (x(i) >= b0 && x(i) <= b1)
            {
                for (int j = 0; j < M; j++)
                {
                    if (y(j) <= s0 && y(j) >= s1)
                    {
                        V(j,i) = V0;
                    }
                    if (y(j) >= s2 || y(j) <= s3)
                    {
                        V(j,i) = V0;
                    }
                }
            }
        }
    }
    if (number_of_slits == 3)
    {
        double s0 = 0.5 + slit_size/2 + slit_size;
        double s1 = 0.5 - slit_size/2 + slit_size;
        double s2 = 0.5 + slit_size/2 - slit_size;
        double s3 = 0.5 - slit_size/2 - slit_size;
        double s4 = 0.5 + 5.0*slit_size/2;
        double s5 = 0.5 - 5.0*slit_size/2;
        for (int i = 0; i < M; i++)
        {
            if (x(i) >= b0 && x(i) <= b1)
            {
                for (int j = 0; j < M; j++)
                {
                    if (y(j) <= s0 && y(j) >= s1)
                    {
                        V(j,i) = V0;
                    }
                    if (y(j) <= s2 && y(j) >= s3)
                    {
                        V(j,i) = V0;
                    }
                    if (y(j) >= s4 || y(j) <= s5)
                    {
                        V(j,i) = V0;
                    }
                }
            }
        }
    }
    V.save("Potential_map.txt",raw_ascii);
}

void Slit_simulation::find_A_and_B()
{
    r = 1i*dt/(2*h*h);
    for (int j = 0; j < M-2; j++)
    {
        for (int i = 0; i < M-2; i++)
        {
            int k = find_k(i,j);
            a(k) = double(1) + double(4)*r + 1i*dt/double(2)*V(i+1,j+1);
            b(k) = double(1) - double(4)*r + 1i*dt/double(2)*V(i+1,j+1);
        }
    }
    A = create_diagonal_matrix(a,-r);
    B = create_diagonal_matrix(b,r);
    /*A.print(cout);
    cout << "\n\n\n" << endl;
    B.print(cout);*/
}

void Slit_simulation::normalise_u()
{
    cx_vec u_conj = conj(u);
    complex<double> norm = sum(cdot(conj(u),u));
    u = u/sqrt(norm);
}

void Slit_simulation::update_U(int slice)
{  
    for (int j = 1; j < M-1; j++)
    {
        for (int i = 1; i < M-1; i++)
        {
            int k = find_k(i-1,j-1);
            U(i,j,slice) = u(k);
        }
    }
}

void Slit_simulation::initial_state(double xc, double yc, double px, double py, double sigmax, double sigmay)
{
    for (int j = 0; j < M-2; j++)
    {
        for (int i = 0; i < M-2; i++)
        {
            int k = find_k(i,j);
            u(k) = exp(-((x(i)-xc)*(x(i)-xc))/(2*sigmax*sigmax) - ((y(j)-yc)*(y(j)-yc))/(2*sigmay*sigmay)
            + 1i*px*((x(i)-xc)) + 1i*py*(y(j)-yc));
        }
    }
    normalise_u();
    print_norm();

    U.zeros(M,M,n);
    update_U(0);
    //u = normalise(u);
}

cx_vec Slit_simulation::compute_b()
{
    cx_vec b_vec(N, fill::zeros);
    for (int i = 0; i < N; i += 3)
    {
        b_vec(i) = (b(i) + double(3)*r)*u(i);
        b_vec(i+1) = (b(i+1) + double(4)*r)*u(i+1);
        b_vec(i+2) = (b(i+2) + double(3)*r)*u(i+2);
    }
    b_vec(0) -= r*u(0);
    b_vec(1) -= r*u(1);
    b_vec(2) -= r*u(2);
    b_vec(N-3) -= r*u(N-3);
    b_vec(N-2) -= r*u(N-2);
    b_vec(N-1) -= r*u(N-1);

    return b_vec;
}

void Slit_simulation::evolve_next_time_step()
{
    for (int i = 0; i < n; i++)
    {
        cx_vec b_vec = B*u;
        //cx_vec b_vec = compute_b();
        u = solve(A,b_vec, solve_opts::likely_sympd);
        update_U(i);
    }
    U.save("Output.dat",raw_ascii);
}