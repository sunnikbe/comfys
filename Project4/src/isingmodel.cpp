
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

Isingmodel::Isingmodel(mat spin_in, double T_in, int seed_in)
{
    spin = spin_in;
    T = T_in;
    L = spin_in.n_rows; // number of particles along one side of the lattice
    E_before = compute_total_energy(spin);

    // generate the RNGs used, the int RNG is for picking new spin states
    // and the real RNG is for accepting the new spin state
    generator.seed(seed_in);
    uniform_int = uniform_int_distribution<int>(0,L-1);
    uniform_dist = uniform_real_distribution<double>(0.0,1.0);

    // compute exponentials
    dE4p = exp(4/T);
    dE4m = exp(-4/T);
    dE8p = exp(8/T);
    dE8m = exp(-8/T);
}

// print out the current spin state
void Isingmodel::print_state()
{
    spin.print(cout);
}

// compute which probability factor used for acceptance
double Isingmodel::compute_prob_factor(int delta_E)
{
    double p_div = 0;
    if (fabs(delta_E) == 8)
    {
        if (delta_E > 0)
        {
            p_div = dE8m;
        }
        else
        {
            p_div = dE8p;
        }
    }
    if (fabs(delta_E) == 4)
    {
        if (delta_E > 0)
        {
            p_div = dE4m;
        }
        else
        {
            p_div = dE4p;
        }
    }
    if (delta_E == 0)
    {
        p_div = 1;
    }
    return p_div;
}

// compute the total energy of the input spin state
int Isingmodel::compute_total_energy(mat spin_)
{
    int E = 0;
    for (int j = 0; j < L; j++)
    {
        for (int i = 0; i < L; i++)
        {
            E += spin_(i,j)*spin_((i+1)%L,j);
            E += spin_(i,j)*spin_(i,(j+1)%L);
        }
    }
    return -E;
}

// compute the total magnetization for the input spin state
int Isingmodel::compute_total_magnetization(mat spin_)
{
    int M = 0;
    for (int j = 0; j < L; j++)
    {
        for (int i = 0; i < L; i++)
        {
            M += spin_(i,j);
        }
    }
    return M;
}

// check if the new spin candidate is accepted or not
void Isingmodel::spin_candidate()
{
    // change one of the spins
    int i = uniform_int(generator);
    int j = uniform_int(generator); 
    // value between 0 and 1 that decides
    // if the spin candidate is accepted
    double r = uniform_dist(generator);
    mat spin_new = spin;
    spin_new(i,j) = -spin(i,j); // new spin candidate

    // compute the difference in total energy delta_E = E_spin_candidate - E_old_spin
    int E_after = compute_total_energy(spin_new);
    int dE = E_after - E_before;
    double p = compute_prob_factor(dE);
    double A = min(1.0,p);

    if (r <= A)
    {
        spin = spin_new;
        E_before = E_after;
    }
}

// Run input amount of MarkovChainMonteCarlo iterations with each MCMC iterations being
// N spin changes, where N is total amount of spins in the system
// writes out the energy per spin and magnetization per spin to file
void Isingmodel::MarkovChainMonteCarlo(string filename, int iterations)
{
    int N = L*L;
    double epsilon = 1.0*E_before/N;
    double m = 1.0*fabs(compute_total_magnetization(spin))/N;
    ofstream ofile(filename, std::ios::out | std::ios::app);
    ofile << setw(5) << epsilon << setw(15) << m << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < iterations; j++)
        {
        spin_candidate();
        epsilon = 1.0*E_before/N;
        m = 1.0*fabs(compute_total_magnetization(spin))/N;
        ofile << setw(5) << epsilon << setw(15) << m << endl;
        }
    }
}