#ifndef MOLECULAR_DYNAMICS_SIMULATION_HPP
#define MOLECULAR_DYNAMICS_SIMULATION_HPP

#include <iostream>
#include <fstream>
#include "ProgressBar.hpp"
#include "ParticleEnsemble.hpp"

class MolecularDynamicsSimulation {
public:
    enum IntegrationMethod {
        Euler,
        EulerCromer,
        SpeedVerlet,
        ThermoEulerCromer,
        ThermoSpeedVerlet,
        NPTEulerCromer,
        NPTSpeedVerlet
    };

private:
    ParticleEnsemble ensemble;
    IntegrationMethod integrationMethod;
    double timeStep;
    int totalSteps;
    double temperature, pressure;
    double taup;
    std::string outputFileName;

public:
    MolecularDynamicsSimulation(int numParticles, double boxLength, double dt, 
                                int numSteps, const std::string& fileName, 
                                IntegrationMethod method, double T, double p, double tau)

    : ensemble(numParticles, boxLength), integrationMethod(method), 
      timeStep(dt), totalSteps(numSteps), temperature(T), pressure(p), 
      taup(tau), outputFileName(fileName) {}


    ~MolecularDynamicsSimulation() {}

    void runEuler(std::ofstream &trajectoryFile);
    void runEulerCromer(std::ofstream &trajectoryFile); 
    void runThermoEulerCromer(std::ofstream &trajectoryFile);   
    void runSpeedVerlet(std::ofstream &trajectoryFile);
    void runThermoSpeedVerlet(std::ofstream &trajectoryFile);
    void runNPTEulerCromer(std::ofstream &trajectoryFile);
    void runNPTSpeedVerlet(std::ofstream &trajectoryFile);

    void run();
};

#endif // MOLECULAR_DYNAMICS_SIMULATION_HPP