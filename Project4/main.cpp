
#include "isingmodel.hpp"
#include <chrono>

using namespace std;
using namespace arma;

int main()
{
    mat spin(20,20, fill::none); // NxN lattice of spins
    double T = 2.4;
    int seed = chrono::system_clock::now().time_since_epoch().count();
    cout << "seed: " << seed << endl;
    cout << endl;

    // fill the spin matrix with 1 and -1 with uniform RNG
    mt19937 generator;
    generator.seed(seed);
    uniform_int_distribution<int> uniform_RNG(1,2);
    spin.imbue( [&]() { return pow(-1,uniform_RNG(generator)); });

    Isingmodel model = Isingmodel(spin,T,seed);
    std::string filename = "Output.txt";
    std::ofstream ofile;
    ofile.open(filename, ofstream::out | ofstream::trunc);
    ofile << setw(-1) << "#" << setw(8) << "epsilon" << setw(11) << "m" << endl; 
    model.print_state();
    model.MarkovChainMonteCarlo(filename,1000);
    ofile.close();

    return 0;
}