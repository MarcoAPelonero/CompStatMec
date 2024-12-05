#include "monteCarloSimulation.hpp"
#include "progressBar.hpp"
#include <iostream>

MonteCarloSimulation::MonteCarloSimulation(int N, ntype T, ntype L, int steps)
    : ensemble(N, L), mcMove(ensemble, T, L), numSteps(steps) {
    // Initialize other members as needed
    ensemble.initializeRandom(L); // Ensure the ensemble is initialized
    ensemble.setInitialEnergy();   // Calculate the initial energy
}

MonteCarloSimulation::~MonteCarloSimulation() {
    if (energyFile.is_open()) {
        energyFile.close();
    }
    if (positionFile.is_open()) {
        positionFile.close();
    }
    if (velocityFile.is_open()) {
        velocityFile.close();
    }
}

void MonteCarloSimulation::initRNG() {
    rng.rseed();
    particleExtractor.rseed();
}

void MonteCarloSimulation::run() {
    initRNG();

    energyFile.open("energies.dat");
    positionFile.open("position.dat");

    if (!energyFile.is_open()) {
        std::cerr << "Error: Could not open energies.dat for writing." << std::endl;
        return;
    }

    if (!positionFile.is_open()) {
        std::cerr << "Error: Could not open position.dat for writing." << std::endl;
        return;
    }

    energyFile << "# Step\tEnergy\n";
    positionFile << "# Step\tPositions\n";

    ProgressBar progressBar(numSteps, 50, "Simulation Progress", "");

    int snapshotInterval = 10; // Save snapshot every 100 steps

    for (int step = 0; step < numSteps; ++step) {
        mcMove.MetropolisStep();
        Energy = ensemble.getEnergy();

        energyFile << step << "\t" << Energy << "\n";

        if (step % snapshotInterval == 0) {
            ensemble.saveSnapshotPosition(positionFile, step);
        }

        progressBar.step();
    }
    mcMove.printAcceptanceRatio();
    energyFile.close();
    positionFile.close();
    progressBar.finish();
}