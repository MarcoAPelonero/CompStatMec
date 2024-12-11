#include <iostream>
#include <cstdlib>
#include <cmath>
#include "monteCarloSimulation.hpp"

int main(int argc, char* argv[]) {
    int numParticles = 200;
    double density = 0.8;   // Default density (in reduced units)
    double temperature = 0.8;
    double delta = 0.5;
    int numSteps = 10000;
    int snapshotInterval = 100; // Save a snapshot every 1000 steps
    int runId = 1;  
    // Optionally parse command-line arguments to customize parameters
    if (argc > 1) numParticles = std::atoi(argv[1]);
    if (argc > 2) density = std::atof(argv[2]);  // Use density as input
    if (argc > 3) temperature = std::atof(argv[3]);
    if (argc > 4) delta = std::atof(argv[4]);
    if (argc > 5) numSteps = std::atoi(argv[5]);
    if (argc > 6) snapshotInterval = std::atoi(argv[6]);
    if (argc > 7) runId = std::atoi(argv[7]);

    // Debugging: Print the parsed values
    std::cout << "numParticles: " << numParticles << std::endl;
    // std::cout << "density: " << density << std::endl;
    // std::cout << "temperature: " << temperature << std::endl;
    // std::cout << "delta: " << delta << std::endl;
    // std::cout << "numSteps: " << numSteps << std::endl;
    // std::cout << "snapshotInterval: " << snapshotInterval << std::endl;
    // std::cout << "runId: " << runId << std::endl;
    // Calculate box length from density and number of particles
    std::cout << "Simulation start" << std::endl;   
    double boxLength = std::pow(numParticles / density, 1.0 / 3.0);
    // Now create the simulation with the calculated box size
    MonteCarloSimulation simulation(numParticles, boxLength, temperature, delta, numSteps, density);
    simulation.runSimulation(snapshotInterval, runId);
    simulation.printAcceptanceRatio();
    std::cout << "Final energy: " << simulation.getEnergy() << std::endl;
    return 0;
}