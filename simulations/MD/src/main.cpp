#include "molecularDynamicsSimulation.hpp"
#include <iostream>
#include <fstream>
#include <cstring>

int main(int argc, char* argv[]) {
    int numParticles = 200; //300
    ntype density = 0.6;
    ntype dt = 0.01;
    int numSteps = 1000;
    ntype temperature = 1;
    std::string outFileName = "ensemble_data.txt";

    if (argc > 1) {
        numParticles = std::stoi(argv[1]);
    }
    if (argc > 2) {
        density = std::stod(argv[2]);
    }
    if (argc > 3) {
        dt = std::stod(argv[3]);
    }
    if (argc > 4) {
        numSteps = std::stoi(argv[4]);
    }
    if (argc > 5) {
        outFileName = argv[5];
    }

    ntype boxSize = std::pow(numParticles / density, 1.0/3.0);
    std::cout << "Box size: " << boxSize << std::endl;
    std::cout << "Out: " << outFileName << std::endl;
    std::ofstream outFile(outFileName);
    if(!outFile.is_open()) {
        std::cerr << "Error opening output file.\n";
        return 1;
    }

    molecularDynamicsSimulation simulation(numParticles, boxSize, dt, numSteps, temperature);
    simulation.setIntegrationMethod(molecularDynamicsSimulation::VELOCITY_VERLET);
    simulation.run(outFile);
    simulation.printAverageEnergy();
    outFile.close();

    return 0;
}
