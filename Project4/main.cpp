
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

int main(int argc, const char* argv[])
{
    if (argc == 1)
    {
        omp_set_num_threads(1); // number of threads you want to use
        int cycles_per_thread = 1e6; // cycles per thread per mcmc cycle
        int seed = chrono::system_clock::now().time_since_epoch().count();

        double T = 1.0; // declare the temperature in units [J/k_B]
        int L = 2; // declare the size of the lattice
        string filename = "set1";

        // create an instance of the Isingmodel class with temperature T [J/k_B]
        // and lattice size L, can add a false statement at the end for
        // randomly generated spin state, if left empty or true the spin state will
        // exist only of +1 spins
        Isingmodel model = Isingmodel(T,L,seed);
        // compute energies and magnetization, set number of threads use with the
        // omp_set_num_thread(N) command where N is number of threads
        model.mcmc(filename,cycles_per_thread);

        // compute energies and magnetization for T = 1.0 and T = 2.4
        // for ordered and unordered starting spin configurations
        string filename_1 = "set1_0_ordered";
        string filename_2 = "set1_0_unordered";
        string filename_3 = "set2_4_ordered";
        string filename_4 = "set2_4_unordered";
        L = 20;
        cycles_per_thread = 3e4;
        Isingmodel model1_0_ordered = Isingmodel(T,L,seed);
        Isingmodel model1_0_unordered = Isingmodel(T,L,seed,false);
        T = 2.4;
        Isingmodel model2_4_ordered = Isingmodel(T,L,seed);
        Isingmodel model2_4_unordered = Isingmodel(T,L,seed,false);
        model1_0_ordered.mcmc(filename_1,cycles_per_thread);
        model1_0_unordered.mcmc(filename_2,cycles_per_thread);
        model2_4_ordered.mcmc(filename_3,cycles_per_thread);
        model2_4_unordered.mcmc(filename_4,cycles_per_thread);

        // compute one million MCMC samples for temperatures T = 1.0 and T = 2.4
        // for ordered starting spin configurations
        string filenamehist1 = "hist_1_0";
        string filenamehist2 = "hist_2_4";
        cycles_per_thread = 1e6;
        Isingmodel model1 = Isingmodel(1.0,L,seed);
        Isingmodel model2 = Isingmodel(2.4,L,seed);
        model1.mcmc(filenamehist1,cycles_per_thread);
        model2.mcmc(filenamehist2,cycles_per_thread);
        
        return 0;
    }

    // Compute energy, magnetization, heat capacity and susceptiblity with the given inputs
    // lattice size L, cycles per thread, burn-in, number of temp steps,
    // number of threads and filename.
    // will use unordered starting spin configurations
    // and will accumulate all the sampled values up to the final cycle
    // with temperatures from 2.1 to 2.4
    if (argc == 7)
    {
        int L = atoi(argv[1]); // lattice size
        int cycles_per_thread = atoi(argv[2]); // cycles per MCMC attempt
        int burn_in = atoi(argv[3]); // burn-in cycles we want to skip
        int number_of_temp_steps = atoi(argv[4]); // number of steps for the temperature
        int number_of_threads = atoi(argv[5]); // amount of threads used
        string filename = argv[6]; // name of file the data is written to
        ofstream ofile;
        // overwrite any data from file with same name
        ofile.open(filename, std::ofstream::out | std::ofstream::trunc);

        // set a random seed from the system clock
        int seed = chrono::system_clock::now().time_since_epoch().count();
        // temperatures used
        vec temps = linspace(2.1,2.4,number_of_temp_steps);
        omp_set_num_threads(number_of_threads);
        // rund n cycles of MCMC for each temperature
        for (double Ts : temps)
        {
            Isingmodel model = Isingmodel(Ts,L,seed,false);
            model.mcmc(filename,cycles_per_thread,burn_in,false);
        }
        ofile.close();

        return 0;
    }
    else
    {
        cout << "need either none or six arguments" << endl;
        exit(1);      
    }
}