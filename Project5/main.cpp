
#include "simulate_slit_experiment.hpp"

using namespace std;
using namespace arma;

int main(int arcg, const char* argv[])
{
    // terminal inputs
    int number_of_slits = atoi(argv[1]); // number og slits used
    double v0 = atof(argv[2]); // strenght of the potential wall
    // set v0 to zero if you don't want any slits or walls
    double sigmay = atof(argv[3]); // input sigma_y (width of the wave packet in y-direction)
    double T = atof(argv[4]); // input total time of simulation

    // set variables used for the simulation and creation of the wavepacket
    double h = 0.005; // positionstep
    int M = 1/h+1; // size of the position plane
    double dt = 2.5e-5; // timestep
    double xc = 0.25; // coordinate of the center of the wave along the x-axis
    double sigmax = 0.05; // width of the wave in x direction
    double px = 200; // momentum in the x direction
    double yc = 0.5; // coordinate of the center of the wave along the y-axis
    double py = 0; // momentum in the y direction, for a wave to hit the slits

    // create an instance of the class with of a given matrix size,
    // time T, timestep dt and position step dx=dy=h
    Slit_simulation slit_sim =  Slit_simulation(M,h,dt,T,v0,number_of_slits);

    // run the class function that creates the potential wall
    //test.compute_potential(v0,thickness,centre,length,slit_size,number_of_slits);
    // set up the initial state for the system
    slit_sim.initial_state(xc,yc,px,py,sigmax,sigmay);
    // evolve the current state for a time T with timestep dt
    slit_sim.evolve_next_time_step();

    return 0;
}
