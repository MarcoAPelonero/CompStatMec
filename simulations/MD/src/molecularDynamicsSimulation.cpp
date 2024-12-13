#include "molecularDynamicsSimulation.hpp"
#include <iomanip>

molecularDynamicsSimulation::molecularDynamicsSimulation(int numParticles, ntype boxSize, ntype dt, int numSteps, ntype temperature)
    : ensemble(numParticles, boxSize), dt(dt), numSteps(numSteps), temperature(temperature), method(EULER) {}

molecularDynamicsSimulation::~molecularDynamicsSimulation() {}

ntype molecularDynamicsSimulation::computeAverageEnergy() {
    ntype totalEnergy = 0.0;
    for (int i = 0; i < ensemble.getNumParticles(); ++i) {
        totalEnergy += ensemble(i).getEnergy();
    }
    return totalEnergy / ensemble.getNumParticles();
}

void molecularDynamicsSimulation::printAverageEnergy() {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Average energy: " << computeAverageEnergy() << std::endl;
}

void molecularDynamicsSimulation::run(std::ofstream &outFile) {
    ensemble.ensembleSnapshot(outFile, true);

    switch (method) {
        case EULER:
            Euler(outFile);
            break;
        case EULER_CROMER:
            EulerCromer(outFile);
            break;
        case LEAP_FROG:
            LeapFrog(outFile);
            break;
        case VELOCITY_VERLET:
            VelocityVerlet(outFile);
            break;
        case THERMO_VELOCITY_VERLET:    
            ThermoVelocityVerlet(outFile);
            break;
    }
}

void molecularDynamicsSimulation::Euler(std::ofstream &outFile) {
    std::cout << "Running Euler" << std::endl;
    ProgressBar bar(numSteps);
    
    for (int i = 0; i < numSteps; ++i) {
        ensemble.updateEnsemble();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepEuler(j, dt);
        }
        bar.update(i);
        ensemble.ensembleSnapshot(outFile);
    }
    bar.finish();
}

void molecularDynamicsSimulation::EulerCromer(std::ofstream &outFile) {
    std::cout << "Running Euler-Cromer" << std::endl;
    ProgressBar bar(numSteps);
    
    for (int i = 0; i < numSteps; ++i) {
        ensemble.updateEnsemble();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepEulerCromer(j, dt);
        }
        bar.update(i);
        ensemble.ensembleSnapshot(outFile);
    }
    bar.finish();
}

void molecularDynamicsSimulation::LeapFrog(std::ofstream &outFile) {
    std::cout << "LeapFrog method not implemented yet." << std::endl;
    ProgressBar bar(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        ensemble.updateEnsemble();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepLeapFrog(j, dt);
        }
        bar.update(i);
        ensemble.ensembleSnapshot(outFile);
    }
    bar.finish();
}

void molecularDynamicsSimulation::VelocityVerlet (std::ofstream &outFile) {
    std::cout << "Running Velocity Verlet" << std::endl;
    int thermoSteps = 10;
    ntype kb = 1.0;
    ProgressBar bar(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        ensemble.updateEnsemble();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepVerlet(j, dt);
        }
        ensemble.updateEnsemblePositions();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepVelocityVerlet(j, dt);
        }
        bar.update(i);
        ensemble.ensembleSnapshot(outFile);
    }
    bar.finish();
}

void molecularDynamicsSimulation::ThermoVelocityVerlet (std::ofstream &outFile) {
    std::cout << "Running Velocity Verlet" << std::endl;
    int thermoSteps = 10;
    ntype kb = 1.0;
    ProgressBar bar(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        ensemble.updateEnsemble();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepVerlet(j, dt);
        }
        ensemble.updateEnsemblePositions();
        for (int j = 0; j < ensemble.getNumParticles(); ++j) {
            ensemble.stepVelocityVerlet(j, dt);
        }
        bar.update(i);
        if (i % thermoSteps == 0) {
            int N = ensemble.getNumParticles();
            ntype totalKinetic = 0.0;

            // Compute the instantaneous kinetic energy
            for (int p = 0; p < N; ++p) {
                totalKinetic += ensemble(p).getKinetic();
            }

            // Compute the current instantaneous temperature (3D case)
            ntype currentTemp = (2.0/3.0) * (totalKinetic / (N * kb));

            // Compute the scaling factor
            ntype scaleFactor = std::sqrt(temperature / currentTemp);

            // Rescale velocities
            for (int p = 0; p < N; ++p) {
                Vector v = ensemble(p).getVelocity();
                v = v * scaleFactor;
                ensemble(p).setVelocity(v);
            }
        }
        // std::cout << std::endl;
        ensemble.ensembleSnapshot(outFile);
    }
    bar.finish();
}

void molecularDynamicsSimulation::setIntegrationMethod(IntegrationMethod method) {
    this->method = method;
}