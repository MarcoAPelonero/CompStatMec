#include "molecularDynamicsSimulation.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
    int numParticles = 1000;
    ntype boxSize = 10.0;
    ntype dt = 0.01;
    int numSteps = 1000;

    if (argc > 1) {
        numParticles = std::stoi(argv[1]);
    }
    if (argc > 2) {
        boxSize = std::stod(argv[2]);
    }
    if (argc > 3) {
        dt = std::stod(argv[3]);
    }
    if (argc > 4) {
        numSteps = std::stoi(argv[4]);
    }

    std::ofstream outFile("ensemble_data.txt");
    if(!outFile.is_open()) {
        std::cerr << "Error opening output file.\n";
        return 1;
    }

    molecularDynamicsSimulation simulation(numParticles, boxSize, dt, numSteps);
    simulation.setIntegrationMethod(molecularDynamicsSimulation::EULER);
    simulation.run(outFile);

    outFile.close();

    return 0;
}
