#include "monteCarloSimulation.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    int N = 2000;           // Number of particles
    double T = 1.0;        // Temperature
    double L = 10;       // Lattice size
    int steps = 1000;       // Number of Monte Carlo steps

    std::cout << "Welcome to the Monte Carlo simulation!" << std::endl;

    if (argc == 5) {
        N = std::stoi(argv[1]);
        T = std::stod(argv[2]);
        L = std::stod(argv[3]);
        steps = std::stoi(argv[4]);
    } else {
        std::cout << "Using default parameters. To customize, provide N, T, L, steps." << std::endl;
    }

    MonteCarloSimulation simulation(N, T, L, steps);
    simulation.run();

    std::cout << "Simulation completed. Energies logged in energies.dat." << std::endl;

    return 0;
}
