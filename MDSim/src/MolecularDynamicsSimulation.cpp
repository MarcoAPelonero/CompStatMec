#include "MolecularDynamicsSimulation.hpp"

void MolecularDynamicsSimulation::run() {

    ProgressBar progressBar(totalSteps);
    // Open trajectory.dat file for writing
    std::ofstream trajectoryFile("trajectory.dat");
    if (!trajectoryFile.is_open()) {
        std::cerr << "Error: could not open trajectory.dat for writing." << std::endl;
        return;
    }

    for (int i = 0; i < totalSteps; ++i) {
        if (i % 50 == 0) {
            ensemble.ensembleThermalize(1.0);
        }
        ensemble.ensembleStepEulerCromer(timeStep);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish();  
    trajectoryFile.close(); 
}