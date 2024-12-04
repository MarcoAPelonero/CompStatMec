#include "monteCarloSimulation.hpp"
#include <iostream>
// ... other includes

MonteCarloSimulation::MonteCarloSimulation(int N, ntype T, ntype L, int steps)
    : ensemble(N), potential(), mcMove(ensemble, T), numSteps(steps), L(L) {
    ensemble.initializeSquareLattice(L);
}

MonteCarloSimulation::~MonteCarloSimulation() {
    if (energyFile.is_open()) {
        energyFile.close();
    }
}

void MonteCarloSimulation::initRNG() {
    rng.rseed();
    particleExtractor.rseed();
}

void MonteCarloSimulation::run() {
    initRNG();

    energyFile.open("energies.dat");
    if (!energyFile.is_open()) {
        std::cerr << "Error: Could not open energies.dat for writing." << std::endl;
        return;
    }

    energyFile << "# Step\tEnergy\n";

    for (int step = 0; step < numSteps; ++step) {
        mcMove.MetropolisStep();
        Energy = ensemble.getEnergy();

        energyFile << step << "\t" << Energy << "\n";
    }

    energyFile.close();
}