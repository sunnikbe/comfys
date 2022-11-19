
#include "isingmodel.hpp"

using namespace std;
using namespace arma;

int main()
{
    omp_set_num_threads(1); // number of threads you want to use
    int cycles_per_thread = 100000; // cycles per thread per mcmc cycle
    
    mat spin_ordered(2,2, fill::ones); // NxN laatice of spins with all spin +1
    mat spin_unordered(20,20, fill::none); // NxN lattice of spins with random spins
    int seed = chrono::system_clock::now().time_since_epoch().count();
    cout << "seed: " << seed << endl;
    cout << endl;

    // fill the spin matrix with 1 and -1 with uniform RNG
    mt19937 generator;
    generator.seed(seed);
    uniform_int_distribution<int> uniform_RNG(1,2);
    spin_unordered.imbue( [&]() { return pow(-1,uniform_RNG(generator)); });

    string filename = "Output.txt";
    ofstream ofile;
    ofile.open(filename, ofstream::out | ofstream::trunc);
    Isingmodel model = Isingmodel(spin_ordered,2.4,seed);
    model.compute_expected_values(filename, cycles_per_thread);

    /*Isingmodel model1_0_ordered = Isingmodel(spin_ordered,1.0,seed);
    Isingmodel model1_0_unordered = Isingmodel(spin_unordered,1.0,seed);
    Isingmodel model2_4_ordered = Isingmodel(spin_ordered,2.4,seed);
    Isingmodel model2_4_unordered = Isingmodel(spin_unordered,2.4,seed);
    string filename1 = "Output_1_0_ordered.txt";
    string filename2 = "Output_1_0_unordered.txt";
    string filename3 = "Output_2_4_ordered.txt";
    string filename4 = "Output_2_4_unordered.txt";
    ofstream ofile1,ofile2,ofile3,ofile4;
    ofile1.open(filename1, ofstream::out | ofstream::trunc);
    ofile2.open(filename2, ofstream::out | ofstream::trunc);
    ofile3.open(filename3, ofstream::out | ofstream::trunc);
    ofile4.open(filename4, ofstream::out | ofstream::trunc);
    //model.print_state();
    vec cycles = {100,500,1000,2000,5000,10000,50000,100000};
    for (int cycle : cycles)
    {
        model1_0_ordered.compute_expected_values(filename1,cycle);
        model1_0_unordered.compute_expected_values(filename2,cycle);
        model2_4_ordered.compute_expected_values(filename3,cycle);
        model2_4_unordered.compute_expected_values(filename4,cycle);

        model1_0_ordered.set_spin_state(spin_ordered);
        model1_0_unordered.set_spin_state(spin_unordered);
        model2_4_ordered.set_spin_state(spin_ordered);
        model2_4_unordered.set_spin_state(spin_unordered);
    }
    ofile1.close();
    ofile2.close();
    ofile3.close();
    ofile4.close();*/

    return 0;
}