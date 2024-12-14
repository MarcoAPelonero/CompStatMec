#include "particleEnsemble.hpp"
#include <cmath>

int main(int argc, char **argv) {
    // Default values
    const char* defaultOutputFile = "output.txt";
    ntype density = 0.6;
    int numParticles = 500;
    int numSteps = 10000;
    ntype dt = 0.001;

    if (argc > 1) {
        defaultOutputFile = argv[1];
    }
    if (argc > 2) {
        density = std::stod(argv[2]);
    }
    if (argc > 3) {
        numParticles = std::stoi(argv[3]);
    }
    if (argc > 4) {
        dt = std::stod(argv[4]);
    }
    if (argc > 5) {
        numSteps = std::stoi(argv[5]);
    }

    std::ofstream outFile(defaultOutputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file " << defaultOutputFile << "\n";
        return 1;
    }

    ntype L = std::cbrt(numParticles / density);
    particleEnsemble ensemble(numParticles, L);

    ensemble.ThermoSpeedVerlet(dt, numSteps, outFile);

    outFile.close();

    return 0;
}