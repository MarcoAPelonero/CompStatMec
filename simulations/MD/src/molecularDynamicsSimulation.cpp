#include "molecularDynamicsSimulation.hpp"

molecularDynamicsSimulation::molecularDynamicsSimulation(int numParticles, ntype boxSize, ntype dt, int numSteps)
    : ensemble(numParticles, boxSize), dt(dt), numSteps(numSteps), method(EULER) {}

molecularDynamicsSimulation::~molecularDynamicsSimulation() {}

void molecularDynamicsSimulation::run(std::ofstream &outFile) {
    ensemble.ensembleSnapshot(outFile, true);

    ProgressBar progressBar(numSteps);

    for (int i = 0; i < numSteps; ++i) {
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            switch (method) {
                case EULER:
                    ensemble.stepEuler(j, dt);
                    break;
                case EULER_CROMER:
                    ensemble.stepEulerCromer(j, dt);
                    break;
            }
        }
        ensemble.updateEnsemble(); 
        ensemble.ensembleSnapshot(outFile, false);
        progressBar.step();
    }
    
    progressBar.finish();
}

void molecularDynamicsSimulation::setIntegrationMethod(IntegrationMethod method) {
    this->method = method;
}