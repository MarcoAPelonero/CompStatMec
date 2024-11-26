#include "sim.hpp"

// The simulation is prepared following a simple cubic lattice as the starting configuration
void simulation::prepare_configuration() {
    int i = 0;
    int nside = retint(cbrt(N_PARTICLES));
    for (int x = 0; x < nside; ++x) {
        for (int y = 0; y < nside; ++y) {
            for (int z = 0; z < nside; ++z) {
                new_particles[i].setPosition(Vector{x*L/N, y*L/N, z*L/N});
                new_particles[i].setSigma(SIGMA);
                new_particles[i].setEpsilon(EPSILON);
                new_particles[i].setRcut(RCUT);
                ++i;
            }
        }
    }
};

void simulation::MonteCarloStep(){

};

void simulation::MonteCarloMetropolisStep(){

};