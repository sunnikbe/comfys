
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

Isingmodel::Isingmodel(double T_in, int L_in, int seed_in, bool ordered_in)
{    
    T = T_in; // temperature of the system
    L = L_in; // number of particles along one side of the lattice
    ordered = ordered_in; // if true the system is ordered, if false it is unordered
    seed = seed_in; // seed for the RNGs

    // create the RNGs
    uniform_int = uniform_int_distribution<int>(0,L-1); // RNG for random spin
    uniform_dist = uniform_real_distribution<double>(0.0,1.0); // RNG for acceptance
    spin_dist = uniform_int_distribution<int>(1,2); // RNG for starting lattice creation
    generator.seed(seed);

    // create the lattice of starting spin configurations
    if (ordered == true)
    {
        spin.ones(L,L);
    }
    else if(ordered == false)
    {
        spin.ones(L,L);
        spin.imbue( [&]() { return pow(-1,spin_dist(generator)); });
    }

    // compute the energy of the system at the start
    E = compute_total_energy(spin); 
    // compute the magnetization of the system at start
    M = compute_total_magnetization(spin);

    // create a container for the ratio of probabilities
    for(int dE = -8; dE <= 8; dE++) prob[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) prob[dE+8] = exp(-dE/T);
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

// run n cycles of Markov Chain Monte Carlo to compute energies and magnetization
// for the Ising model, each cycle is written down into a matrix.
// uses openMP for parallelizing the code
void Isingmodel::mcmc(string filename, int cycles, int burn_in, bool mat_form)
{
    #pragma omp parallel
    {    
        int N = L*L; // number of particles in the system
        int num_threads = omp_get_num_threads();  
        
        // create a matrix for the energies and magnetization
        // would have moved it outside the parallel, but couldn't get it to work
        // and had more important issues
        E_mat.zeros(cycles,num_threads);
        EE_mat.zeros(cycles,num_threads);
        M_mat.zeros(cycles,num_threads);
        MM_mat.zeros(cycles,num_threads);

        int E_thread = E; // private energy variable
        int M_thread = M; // private magnetization variable
        int dE = 0;
        mat spin_thread = spin;
        // if multiple threads used, this ensures that they have different seeds
        int thread_num = omp_get_thread_num();
        generator.seed(seed + thread_num);
        // run MCMC for number of cycles
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

                //dcompute the new delta E
                dE = 2*spin_thread(i,j) *
                    (spin_thread((i+1)%L,j) + spin_thread(i,(j+1)%L) +
                    spin_thread(((i-1)%L+L)%L,j) + spin_thread(i,((j-1)%L+L)%L));
                // check for acceptance
                if (r <= prob[dE+8])
                {
                    spin_thread(i,j) *= -1;
                    E_thread += dE;
                    M_thread += 2*spin_thread(i,j);
                }
                // add the energy and magnetization for current MCMC cycle
                E_mat(l,thread_num) += E_thread;
                EE_mat(l,thread_num) += E_thread*E_thread;
                M_mat(l,thread_num) += fabs(M_thread);
                MM_mat(l,thread_num) += M_thread*M_thread;
            }
        }

    }
    // if the file should be written as the matrices
    if (mat_form == 1)
    {
        E_mat.save("E_mat_"+filename+".txt",raw_ascii);
        EE_mat.save("EE_mat_"+filename+".txt",raw_ascii);
        M_mat.save("M_mat_"+filename+".txt",raw_ascii);
        MM_mat.save("MM_mat_"+filename+".txt",raw_ascii);
    }

    // or rather compute the energy, magnetization, heat capacity and susceptiblity
    // accumulated for all the cycles
    int N = L*L;
    int n_threads = omp_get_max_threads();
    double E = 0;
    double EE = 0;
    double M = 0;
    double MM = 0;
    // can remove the first few instances of samples with burn-in
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

    // compute energies and magnetization per spin
    double eps = E/N;
    double m = M/N;

    // compute heat capacity and susceptiblity
    double Evar = (EE - E*E)/N;
    double Mvar = (MM - M*M)/N;
    double CV = Evar/T/T;
    double chi = Mvar/T;
    // write to file
    if (mat_form == 0)
    {
        ofstream ofile(filename, std::ios::out | std::ios::app);
        ofile << setw(6) << T << 
        setw(14) << eps << 
        setw(14) << m << 
        setw(14) << CV << 
        setw(14) << chi << endl;
    }
    // print to terminal values computed
    cout << L << "x" << L << " lattice, " << "T = " << T << " J/k_B" << endl;
    cout << "E = " << E << "  E^2 = " << EE << endl;
    cout << "|M| = " << M << "  M^2 = " << MM << endl;
    cout << "Epsilon = " << eps << endl;
    cout << "m = " << m << endl;
    cout << "C_V = " << CV << endl;
    cout << "chi = " << chi << endl;
    cout << endl;

}