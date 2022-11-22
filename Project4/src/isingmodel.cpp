
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

Isingmodel::Isingmodel(double T_in, int L_in, int seed_in, bool ordered_in)
{    
    T = T_in;
    L = L_in; // number of particles along one side of the lattice
    ordered = ordered_in;
    seed = seed_in;

    uniform_int = uniform_int_distribution<int>(0,L-1);
    uniform_dist = uniform_real_distribution<double>(0.0,1.0);
    spin_dist = uniform_int_distribution<int>(1,2);
    generator.seed(seed);

    //spin = spin_in;
    if (ordered == true)
    {
        spin.ones(L,L);
    }
    else if(ordered == false)
    {
        spin.ones(L,L);
        spin.imbue( [&]() { return pow(-1,spin_dist(generator)); });
    }

    E = compute_total_energy(spin);
    M = compute_total_magnetization(spin);

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

// run n cycles of Markov Chain Monte Carlo to compute energies and magnetization
// for the Ising model, each cycle is written down into a matrix.
// uses openMP for parallelizing the code
void Isingmodel::mcmc(string filename, int cycles, int burn_in, bool mat_form)
{
    #pragma omp parallel
    {    
        int N = L*L;
        int num_threads = omp_get_num_threads();  
        
        E_mat.zeros(cycles,num_threads);
        EE_mat.zeros(cycles,num_threads);
        M_mat.zeros(cycles,num_threads);
        MM_mat.zeros(cycles,num_threads);

        int E_thread = E;
        int M_thread = M;
        int dE = 0;
        mat spin_thread = spin;
        int thread_num = omp_get_thread_num();
        generator.seed(seed + thread_num);
        for (int l = 0; l < cycles; l++) 
        {
            for (int k = 0; k < N; k++)
            {
                // change one of the spins
                int i = uniform_int(generator);
                int j = uniform_int(generator); 
                // value between 0 and 1 that decides
                // if the spin candidate is accepted
                double r = uniform_dist(generator);

                //dE = compute_delta_E(i,j,spin_thread);
                dE = 2*spin_thread(i,j) *
                    (spin_thread((i+1)%L,j) + spin_thread(i,(j+1)%L) +
                    spin_thread(((i-1)%L+L)%L,j) + spin_thread(i,((j-1)%L+L)%L));
                if (r <= prob[dE+8])
                {
                    spin_thread(i,j) *= -1;
                    E_thread += dE;
                    M_thread += 2*spin_thread(i,j);
                }
                E_mat(l,thread_num) += E_thread;
                EE_mat(l,thread_num) += E_thread*E_thread;
                M_mat(l,thread_num) += fabs(M_thread);
                MM_mat(l,thread_num) += M_thread*M_thread;
            }
        }

    }
    if (mat_form == 1)
    {
        E_mat.save("E_mat_"+filename+".txt",raw_ascii);
        EE_mat.save("EE_mat_"+filename+".txt",raw_ascii);
        M_mat.save("M_mat_"+filename+".txt",raw_ascii);
        MM_mat.save("MM_mat_"+filename+".txt",raw_ascii);
    }
    int N = L*L;
    int n_threads = omp_get_max_threads();
    double E = 0;
    double EE = 0;
    double M = 0;
    double MM = 0;
    for (int i = burn_in; i < cycles; i++)
    {
        E += sum(E_mat.row(i))/n_threads;
        EE += sum(EE_mat.row(i))/n_threads;
        M += sum(M_mat.row(i))/n_threads;
        MM += sum(MM_mat.row(i))/n_threads;
    }
    int ns = cycles-burn_in;
    E = 1.0*E/ns/N;
    EE = 1.0*EE/ns/N;
    M = 1.0*M/ns/N;
    MM = 1.0*MM/ns/N;

    double eps = E/N;
    double m = M/N;

    double Evar = (EE - E*E)/N;
    double Mvar = (MM - M*M)/N;
    double CV = Evar/T/T;
    double xi = Mvar/T;
    if (mat_form == 0)
    {
        ofstream ofile(filename, std::ios::out | std::ios::app);
        ofile << setw(6) << T << 
        setw(14) << eps << 
        setw(14) << m << 
        setw(14) << CV << 
        setw(14) << xi << endl;
    }
    cout << L << "x" << L << " lattice, " << "T = " << T << " J/k_B" << endl;
    cout << "E = " << E << "  E^2 = " << EE << endl;
    cout << "|M| = " << M << "  M^2 = " << MM << endl;
    cout << "Epsilon = " << eps << endl;
    cout << "m = " << m << endl;
    cout << "C_V = " << CV << endl;
    cout << "chi = " << xi << endl;
    cout << endl;

}