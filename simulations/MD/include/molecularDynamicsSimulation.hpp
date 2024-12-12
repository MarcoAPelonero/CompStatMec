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
            VELOCITY_VERLET 
            };

    private:
        particleEnsemble ensemble;
        ntype dt;
        int numSteps;
        IntegrationMethod method;

    public:
        molecularDynamicsSimulation(int numParticles, ntype boxSize, ntype dt, int numSteps);
        ~molecularDynamicsSimulation();

        ntype computeAverageEnergy();
        void printAverageEnergy();

        void run(std::ofstream &outFile);
        void setIntegrationMethod(IntegrationMethod method);
};

#endif // MOLECULARDYNAMICSSIMULATION_HPP
