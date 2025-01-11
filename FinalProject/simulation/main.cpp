#include "MolecularDynamicsSimulation.hpp"
#include <iostream>
#include <cstdlib>
#include <filesystem>

MolecularDynamicsSimulation::IntegrationMethod parseMethod(const std::string& methodStr) {
    if (methodStr == "Euler") return MolecularDynamicsSimulation::IntegrationMethod::Euler;
    if (methodStr == "EulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::EulerCromer;
    if (methodStr == "SpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::SpeedVerlet;
    if (methodStr == "ThermoEulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::ThermoEulerCromer;
    if (methodStr == "ThermoSpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::ThermoSpeedVerlet;
    if (methodStr == "NPTEulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::NPTEulerCromer;
    if (methodStr == "NPTSpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::NPTSpeedVerlet;
    if (methodStr == "RDFEulerCromer") return MolecularDynamicsSimulation::IntegrationMethod::RDFEulerCromer;
    if (methodStr == "RDFSpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::RDFSpeedVerlet;
    if (methodStr == "MLSpeedVerlet") return MolecularDynamicsSimulation::IntegrationMethod::MLSpeedVerlet;
    throw std::invalid_argument("Invalid integration method: " + methodStr);
}

int main(int argc, char* argv[]) {
    // Default parameters
    int N = 300;
    double rho = 0.6;
    double T = 1.0;
    double p = 1;
    double dt = 0.005;
    int numSteps = 10000;

    std::string fileName = "trajectory.dat";
    std::string rdfName = "radialDistribution.dat";
    std::string mlName = "mlComparison.dat";

    double taup = 1;
    std::cout << "Current Working Directory: " << std::filesystem::current_path() << std::endl;
    MolecularDynamicsSimulation::IntegrationMethod method = MolecularDynamicsSimulation::IntegrationMethod::RDFSpeedVerlet;

    // Override with command-line arguments if provided
    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) rho = std::atof(argv[2]);
    if (argc > 3) dt = std::atof(argv[3]);
    if (argc > 4) numSteps = std::atoi(argv[4]);
    if (argc > 5) fileName = argv[5];
    if (argc > 6) rdfName = argv[6];
    if (argc > 7) {
        try {
            method = parseMethod(argv[7]);
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }
    if (argc > 8) T = std::atof(argv[8]);
    if (argc > 9) p = std::atof(argv[9]);
    if (argc > 10) taup = std::atof(argv[10]);
    if (argc > 11) mlName = argv[11];

    double boxLenght = std::pow(N / rho, 1.0 / 3.0);
    
    MolecularDynamicsSimulation simulation(N, boxLenght, dt, numSteps, fileName, method, T, p, taup, rdfName, mlName);

    simulation.run();

    return 0;
}
