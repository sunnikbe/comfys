
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

Isingmodel::Isingmodel(mat spin_in, double T_in, int seed_in)
{
    spin = spin_in;
    T = T_in;
    L = spin_in.n_rows; // number of particles along one side of the lattice
    E = compute_total_energy(spin);
    M = compute_total_magnetization(spin);

    uniform_int = uniform_int_distribution<int>(0,L-1);
    uniform_dist = uniform_real_distribution<double>(0.0,1.0);
    seed = chrono::system_clock::now().time_since_epoch().count();

    for(int dE = -8; dE <= 8; dE++) prob[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) prob[dE+8] = exp(-dE/T);
}

// print out the current spin state
void Isingmodel::print_state()
{
    spin.print(cout);
}

void Isingmodel::set_spin_state(mat spin_in)
{
    spin = spin_in;
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

int Isingmodel::compute_delta_E(int i, int j, mat spin_in)
{
    return 2*spin_in(i,j) *
    (spin_in((i+1)%L,j) + spin_in(i,(j+1)%L) +
    spin_in(((i-1)%L+L)%L,j) + spin_in(i,((j-1)%L+L)%L));
}

// Run input amount of MarkovChainMonteCarlo iterations with each MCMC iterations being
// N spin changes, where N is total amount of spins in the system
// writes out the energy per spin and magnetization per spin to file
void Isingmodel::MCMC(int index, int thread_id,int E_thread, int M_thread, int dE,
mat spin_thread)
{
    int N = L*L;
    for (int k = 0; k < N; k++)
    {
        // change one of the spins
        int i = uniform_int(generator);
        int j = uniform_int(generator); 
        // value between 0 and 1 that decides
        // if the spin candidate is accepted
        double r = uniform_dist(generator);

        dE = compute_delta_E(i,j,spin_thread);
        if (fabs(dE) != 8 && fabs(dE) != 4 && dE != 0)
        {
            #pragma omp master
            {
                cout << dE << "  " << i << "  " << j << endl;
            }
        }
        if (r <= prob[dE+8])
        {
            spin_thread(i,j) *= -1;
            E_thread += dE;
            M_thread += 2*spin_thread(i,j);
        }
        E_mat(index,thread_id) += 1.0*E_thread/N;
        EE_mat(index,thread_id) += 1.0*E_thread*E_thread/N;
        M_mat(index,thread_id) += 1.0*fabs(M_thread)/N;
        MM_mat(index,thread_id) += 1.0*M_thread*M_thread/N;
    }
}

void Isingmodel::compute_expected_values(string filename, int cycles)
{
    #pragma omp parallel private(prob,generator)
    {    
        int N = L*L;
        int num_threads = omp_get_num_threads();  
        E_mat.zeros(cycles,num_threads);
        EE_mat.zeros(cycles,num_threads);
        M_mat.zeros(cycles,num_threads);
        MM_mat.zeros(cycles,num_threads);

        int E_thread,M_thread;
        int dE = 0;
        mat spin_thread;
        #pragma omp critical
        {
            E_thread = E;
            M_thread = M;
            spin_thread = spin;
        }
        int thread_num = omp_get_thread_num();
        generator.seed(seed + thread_num);
        for (int l = 0; l < cycles; l++) 
        {
            MCMC(l,thread_num,E_thread,M_thread,dE,spin_thread);
        }

    }
    E_mat.save("E_mat.txt",raw_ascii);
    EE_mat.save("EE_mat.txt",raw_ascii);
    M_mat.save("M_mat.txt",raw_ascii);
    MM_mat.save("MM_mat.txt",raw_ascii);
    //E_mat.print(cout);
    cout << endl;
    int N_part = L*L;
    int n_threads = omp_get_max_threads();
    //cout << n_threads << endl;
    cout << sum(E_mat,0)/cycles/N_part << endl;
    cout << accu(E_mat)/cycles/N_part/n_threads << endl;
    cout << endl;
    cout << sum(M_mat,0)/cycles/N_part << endl;
    cout << accu(M_mat)/cycles/N_part/n_threads << endl;

}