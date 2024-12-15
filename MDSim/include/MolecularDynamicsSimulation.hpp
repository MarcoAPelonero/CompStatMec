#ifndef MOLECULAR_DYNAMICS_SIMULATION_HPP
#define MOLECULAR_DYNAMICS_SIMULATION_HPP

#include "progressBar.hpp"
#include "particleEnsemble.hpp"

class MolecularDynamicsSimulation {
    private:
        ParticleEnsemble ensemble;
        double timeStep;
        int totalSteps;
    public:
        MolecularDynamicsSimulation(int numParticles, double boxLength, double dt, int numSteps)
            : ensemble(numParticles, boxLength), timeStep(dt), totalSteps(numSteps) {}
        ~MolecularDynamicsSimulation() {}

        void run();
};

#endif // MOLECULAR_DYNAMICS_SIMULATION_HPP