#include "MolecularDynamicsSimulation.hpp"
#include <iostream>
#include <cstdlib>

MolecularDynamicsSimulation::IntegrationMethod parseMethod(const std::string& methodStr) {
    if (methodStr == "Euler") return MolecularDynamicsSimulation::IntegrationMethod::Euler;
    if (methodStr == "EulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::EulerCromer;
    if (methodStr == "SpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::SpeedVerlet;
    if (methodStr == "ThermoEulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::ThermoEulerCromer;
    if (methodStr == "ThermoSpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::ThermoSpeedVerlet;
    if (methodStr == "NPTEulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::NPTEulerCromer;
    if (methodStr == "NPTSpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::NPTSpeedVerlet;
    throw std::invalid_argument("Invalid integration method: " + methodStr);
}

int main(int argc, char* argv[]) {
    // Default parameters
    int N = 600;
    double rho = 0.6;
    double T = 1.0;
    double p = 1.011;
    double dt = 0.003;
    int numSteps = 5000;
    std::string fileName = "trajectory.dat";
    double taup = 1;

    MolecularDynamicsSimulation::IntegrationMethod method = MolecularDynamicsSimulation::IntegrationMethod::ThermoSpeedVerlet;

    // Override with command-line arguments if provided
    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) rho = std::atof(argv[2]);
    if (argc > 3) dt = std::atof(argv[3]);
    if (argc > 4) numSteps = std::atoi(argv[4]);
    if (argc > 5) fileName = argv[5];
    if (argc > 6) {
        try {
            method = parseMethod(argv[6]);
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }
    if (argc > 7) T = std::atof(argv[7]);
    if (argc > 8) p = std::atof(argv[8]);
    if (argc > 9) taup = std::atof(argv[9]);

    std::cout << "Running simulation with: "<< 
        "\nN = " << N << 
        "\nrho = " << rho << 
        "\nT = " << T << 
        "\ndt = " << dt << 
        "\nnumSteps = " << numSteps << 
        "\nfileName = " << fileName << std::endl;

    double boxLenght = std::pow(N / rho, 1.0 / 3.0);
    
    MolecularDynamicsSimulation simulation(N, boxLenght, dt, numSteps, fileName, method, T, p, taup);

    simulation.run();

    return 0;
}
