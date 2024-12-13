#ifndef MOLECULARDYNAMICSSIMULATION_HPP
#define MOLECULARDYNAMICSSIMULATION_HPP

#include "particleEnsemble.hpp"
#include "progressBar.hpp"

class molecularDynamicsSimulation {
    public:
        enum IntegrationMethod { 
            EULER, 
            EULER_CROMER,
            LEAP_FROG,
            VELOCITY_VERLET,
            THERMO_VELOCITY_VERLET
            };

    private:
        particleEnsemble ensemble;
        ntype dt;
        int numSteps;
        ntype temperature;
        IntegrationMethod method;

    public:
        molecularDynamicsSimulation(int numParticles, ntype boxSize, ntype dt, int numSteps, ntype temperature);
        ~molecularDynamicsSimulation();

        ntype computeAverageEnergy();
        void printAverageEnergy();

        void run(std::ofstream &outFile);
        void setIntegrationMethod(IntegrationMethod method);

        void Euler(std::ofstream &outFile);
        void EulerCromer(std::ofstream &outFile); 
        void LeapFrog(std::ofstream &outFile);
        void VelocityVerlet(std::ofstream &outFile);
        void ThermoVelocityVerlet(std::ofstream &outFile);
};

#endif // MOLECULARDYNAMICSSIMULATION_HPP
