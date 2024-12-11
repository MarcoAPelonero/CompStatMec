#include "monteCarloSimulation.hpp"
#include "progressBar.hpp"  // Progress bar header
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream> // For stringstream

MonteCarloSimulation::MonteCarloSimulation(int numParticles, double boxLength, double T, double delta, int numSteps, ntype dens)
    : ensemble(numParticles, boxLength), mcMove(ensemble, T, delta),
      temperature(T), displacement(delta), steps(numSteps), density(dens) {
    ensemble.initializeEnsemble();
    rng.rseed();
    particleExtractor.setDistribution(0, numParticles - 1);
    particleExtractor.rseed();
}

MonteCarloSimulation::~MonteCarloSimulation() {}

void MonteCarloSimulation::runSimulation(int snapshotInterval, int runId) {

    // Generate unique filenames using runId
    std::stringstream ss_energy;
    ss_energy << "data/energy_" << runId << "_" << temperature << "_" << density << ".txt";
    std::ofstream energyFile(ss_energy.str());

    std::stringstream ss_positions;
    ss_positions << "data/positions_" << runId << "_" << temperature << "_" << density << ".txt";
    std::ofstream snapshotFile(ss_positions.str());

    // Check if files are successfully opened
    if (!snapshotFile.is_open() || !energyFile.is_open()) {
        std::cerr << "Error opening output files for run " << runId << "!" << std::endl;
        return;
    }

    ProgressBar progressBar(steps);
    // Loop over steps and perform the simulation
    for (int i = 1; i <= steps; ++i) {
        mcMove.MetropolisStep();

        // Update the progress bar
        progressBar.step();

        // Save snapshot and energy at the specified intervals
        if (snapshotInterval > 0 && i % snapshotInterval == 0) {
            ensemble.saveSnapshotPosition(snapshotFile, i);
            
            // Get energy and write to file
            double energy = ensemble.getEnergy();
            energyFile << i << " " << energy << std::endl;
        }
    }

    progressBar.finish();
    // Close the output files
    snapshotFile.close();
    energyFile.close();
}

void MonteCarloSimulation::printAcceptanceRatio() const {
    mcMove.printAcceptanceRatio();
}

double MonteCarloSimulation::getAcceptanceRatio() const {
    return mcMove.getAcceptanceRatio();
}

double MonteCarloSimulation::getEnergy() const {
    return ensemble.getEnergy();
}
