
#include "simulate_slit_experiment.hpp"

using namespace std;
using namespace arma;

int main(int arcg, const char* argv[])
{
    int M = atoi(argv[1]);
    int number_of_slits = atoi(argv[2]);

    double h = 0.005;
    double dt = 2.5e-5;
    double T = 0.002;
    double xc = 0.25;
    double sigmax = 0.05;
    double px = 200;
    double yc = 0.5;
    double sigmay = 0.2;
    double py = 0;

    Slit_simulation test =  Slit_simulation(M,h,dt,T);
    test.compute_potential(1e10,0.02,0.5,0.05,0.05,number_of_slits);
    test.initial_state(xc,yc,px,py,sigmax,sigmay);

    //test.print_u();

    test.evolve_next_time_step();
    //test.print_U();
    //test.print_norm();
    return 0;
}