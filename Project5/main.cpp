
#include "simulate_slit_experiment.hpp"

using namespace std;
using namespace arma;

int main()
{
    int M = 5;

    double h = 0.005;
    double dt = 2.5e-5;
    double T = 0.002;
    double xc = 0.25;
    double sigmax = 0.05;
    double px = 200;
    double yc = 0.5;
    double sigmay = 0.2;
    double py = 0;

    Slit_simulation test =  Slit_simulation(M,h,dt);
    test.compute_potential(1e10,0.02,0.5,0.05,0.05,2);
    test.find_A_and_B();
    test.initial_state(xc,yc,px,py,sigmax,sigmay);
    test.evolve_next_time_step();

    return 0;
}