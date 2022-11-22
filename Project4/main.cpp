
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

int main(int argc, const char* argv[])
{
    if (argc == 1)
    {
        omp_set_num_threads(1); // number of threads you want to use
        int cycles_per_thread = 1e5; // cycles per thread per mcmc cycle
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
        
        return 0;
    }

    if (argc == 6)
    {
        int L = atoi(argv[1]); // lattice size
        int cycles_per_thread = atoi(argv[2]); // cycles per MCMC attempt
        int number_of_temp_steps = atoi(argv[3]); // number of steps for the temperature
        int number_of_threads = atoi(argv[4]); // amount of threads used
        string filename = argv[5]; // name of file the data is written to
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
            model.mcmc(filename,cycles_per_thread,100,false);
        }
        ofile.close();

        return 0;
    }
}