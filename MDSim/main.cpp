#include "MolecularDynamicsSimulation.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
    // Default parameters
    int N = 300;
    double rho = 0.6;
    double T = 1.0;
    double dt = 0.001;
    int numSteps = 20000;

    // Override with command-line arguments if provided
    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) rho = std::atof(argv[2]);
    if (argc > 3) T = std::atof(argv[3]);
    if (argc > 4) dt = std::atof(argv[4]);
    if (argc > 5) numSteps = std::atoi(argv[5]);

    double boxLength = std::pow(N / rho, 1.0 / 3.0);

    MolecularDynamicsSimulation sim(N, boxLength, dt, numSteps);
    sim.run();
    return 0;
}