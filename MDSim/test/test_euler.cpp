#include <iostream>
#include "particleEnsemble.hpp" // Adjust path as needed
#include <cmath>
#include <random>

// A simple test to verify that the Euler integration updates positions and velocities.

int main() {
    // Parameters for the test
    int numParticles = 2;          // A small number of particles
    double boxLength = 10.0;       // Size of the simulation box
    double dt = 0.001;             // Time step
    int numSteps = 1;           // Number of integration steps

    // Create a ParticleEnsemble
    ParticleEnsemble ensemble(numParticles, boxLength);

    // Print initial state
    std::cout << "Initial State:" << std::endl;
    ensemble.showParticles();

    ensemble.computeForces();
    ensemble.computePotentials();
    ensemble.computeKinetics();
    ensemble.computeEnergies();
    ensemble.computeEnsembleEnergies();

    double initialEnergy = ensemble.getTotalEnergy();
    std::cout << "Initial Total Energy: " << initialEnergy << std::endl;

    // Run integration steps
    for (int step = 0; step < numSteps; ++step) {
        ensemble.ensembleStepEuler(dt);

        if (step % 100 == 0) {
            double totalEnergy = ensemble.getTotalEnergy();
            double kineticEnergy = ensemble.getTotalKineticEnergy();
            double potentialEnergy = ensemble.getTotalPotentialEnergy();

            std::cout << "Step: " << step 
                      << "  Kinetic Energy: " << kineticEnergy
                      << "  Potential Energy: " << potentialEnergy
                      << "  Total Energy: " << totalEnergy << std::endl;

            // Optionally show a few particles to check position updates
            // ensemble.showParticles();
        }
    }

    // After simulation
    double finalEnergy = ensemble.getTotalEnergy();
    std::cout << "Final Total Energy: " << finalEnergy << std::endl;
    double energyChange = (finalEnergy - initialEnergy)/initialEnergy * 100.0;
    std::cout << "Percentage Energy Change: " << energyChange << "%" << std::endl;

    return 0;
}
