
#include "simulate_slit_experiment.hpp"

using namespace std;
using namespace arma;

Slit_simulation::Slit_simulation(int M_in, double h_in, double dt_in, double T_in,
                                 double v0_in, int number_of_slits_in)
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

    thickness = 0.02;
    centre = 0.5;
    length = 0.05;
    slit_size = 0.05;
    v0 = v0_in;
    number_of_slits = number_of_slits_in;
    V = mat(M,M,fill::zeros);
    if (v0 != 0)
    {
        compute_potential();
    }
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
    cout << scientific << cdot(u,u) << endl;
}

sp_cx_mat Slit_simulation::create_diagonal_matrix(cx_vec d, complex<double> r)
{
    cx_vec r_vec(N-1, fill::value(r));
    for (int i = M-3; i < N-1; i += M-2)
    {
        r_vec(i) = 0;
    }
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
    return i + j*(M-2);
}

void Slit_simulation::compute_potential()
{
    double b0 = centre - thickness;
    double b1 = centre + thickness;
    int n_slit0 = (0.5 - (2*number_of_slits - 1)/2*slit_size)/1.0*M;
    int n_slit1 = (0.5 + (2*number_of_slits - 1)/2*slit_size)/1.0*M;
    int n_dslit = slit_size/h;
    vec v_i(M);
    for (int i = b0/1.0*M; i < b1/1.0*M; i++)
    {
        v_i.fill(v0);
        for (int j = 0; j < number_of_slits; j++)
        {
            int n0 = n_slit0 + 2*j*n_dslit;
            int n1 = n_slit0 + n_dslit + 2*j*n_dslit - 1;
            v_i(span(n0,n1)).zeros();
        }
        V.col(i) = v_i;
        
        /*V.col(i).fill(1.0);
        for (int j = 0; j < number_of_slits; j++)
        {
            V.col(i)(span(n_slit0+2*n_dslit,n_slit1+2*n_dslit)).fill(0);
            for (int k = n_slit0+2*n_dslit; k < n_slit1+2*n_dslit; k++)
            {
                V.col(i)
            }
        }*/
    }

    /*V = mat(M,M, fill::zeros);
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
                        V(i,j) = V0;
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
                    //cout << i << endl;
                    if (y(j) <= s0 && y(j) >= s1)
                    {
                        V(i,j) = V0;
                    }
                    if (y(j) >= s2 || y(j) <= s3)
                    {
                        V(i,j) = V0;
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
                        V(i,j) = V0;
                    }
                    if (y(j) <= s2 && y(j) >= s3)
                    {
                        V(i,j) = V0;
                    }
                    if (y(j) >= s4 || y(j) <= s5)
                    {
                        V(i,j) = V0;
                    }
                }
            }
        }
    }*/
    V.save("Potential_map.txt",raw_ascii);
}

void Slit_simulation::find_A_and_B()
{
    r = 1.0i*dt/double(2)/h/h;
    for (int j = 0; j < M-2; j++)
    {
        for (int i = 0; i < M-2; i++)
        {
            int k = find_k(i,j);
            //cout << k << "  " << i << "  " << j << endl;
            a(k) = double(1) + double(4)*r + 1i*dt/double(2)*V(j+1,i+1);
            b(k) = double(1) - double(4)*r - 1i*dt/double(2)*V(j+1,i+1);
        }
    }
    A = create_diagonal_matrix(a,-r);
    B = create_diagonal_matrix(b,r);
}

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
    complex<double> norm = cdot(u,u);
    u = u/sqrt(norm);
    U.zeros(M,M,n+1);
    update_U(0);
}

void Slit_simulation::evolve_next_time_step()
{
    find_A_and_B();
    superlu_opts opts;
    opts.symmetric = true;
    for (int i = 0; i <= n; i++)
    {
        cx_vec b_vec = B*u;
        spsolve(u,A,b_vec, "superlu", opts);
        update_U(i);
        //print_norm();
    }
    U.save("Output.bin",arma_binary);
}