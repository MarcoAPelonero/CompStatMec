#include "particleEnsemble.hpp"
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <limits>
#include <fstream>

particleEnsemble::particleEnsemble(int numParticles, ntype boxSize)
    : numParticles(numParticles), boxSize(boxSize), potential(boxSize) {
    
    particles.reserve(numParticles);

    int particlesPerSide = static_cast<int>(std::ceil(std::cbrt(static_cast<ntype>(numParticles))));
    ntype spacing = boxSize / particlesPerSide;

    // Seed the random engine with a deterministic seed for reproducibility
    // or use std::random_device for a non-deterministic seed if desired.
    std::default_random_engine generator(42);
    std::normal_distribution<ntype> distribution(0.0, 1.0); 

    int count = 0;
    for (int i = 0; i < particlesPerSide && count < numParticles; ++i) {
        for (int j = 0; j < particlesPerSide && count < numParticles; ++j) {
            for (int k = 0; k < particlesPerSide && count < numParticles; ++k) {
                Particle p(Vector{i * spacing, j * spacing, k * spacing});

                ntype vx = distribution(generator);
                ntype vy = distribution(generator);
                ntype vz = distribution(generator);
                p.setVelocity(Vector{vx, vy, vz});

                p.store();
                particles.push_back(p);
                ++count;
            }
        }
    }

    // Compute minimum initial distance between particles (for diagnostic)
    ntype minDistance = std::numeric_limits<ntype>::max();
    for (size_t i = 0; i < particles.size(); ++i) {
        const Vector &pos_i = particles[i].getPosition();
        for (size_t j = i + 1; j < particles.size(); ++j) {
            ntype distance = (pos_i - particles[j].getPosition()).modulus();
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
    }
    std::cout << "Minimum distance between particles: " << minDistance << std::endl;

    ensembleEnergy = 0.0;
    ensembleStep();
}

particleEnsemble::~particleEnsemble() {
    particles.clear();
}

Particle& particleEnsemble::operator()(int i) {
    return particles[i];
}

int particleEnsemble::getNumParticles() {
    return numParticles;
}

ntype particleEnsemble::getEnsembleEnergy() {
    return ensembleEnergy;
}

void particleEnsemble::show() {
    for (int i = 0; i < numParticles; ++i) {
        std::cout << "Particle " << i+1 << ":\n";
        particles[i].show();
    }
}

inline void particleEnsemble::applyPeriodicBoundary(Vector &pos) {
    for (int i = 0; i < dim; ++i) {
        pos(i) = std::fmod(pos(i), boxSize);
        if (pos(i) < 0) {
            pos(i) += boxSize;
        }
    }
}

void particleEnsemble::ensembleSnapshot(std::ofstream &outFile, bool writeHeader) {
    if (writeHeader) {
        outFile << "x y z vx vy vz energy\n";
    }
    for (int i = 0; i < numParticles; ++i) {
        const Vector &pos = particles[i].getPosition();
        const Vector &vel = particles[i].getVelocity();
        ntype E = particles[i].getPotential();
        outFile << pos(0) << " " << pos(1) << " " << pos(2) << " "
                << vel(0) << " " << vel(1) << " " << vel(2) << " "
                << E << "\n";
    }
    outFile << "\n";
}

void particleEnsemble::particleEnergyStep(int particleIndex){
    const Vector position = particles[particleIndex].getPosition();
    ntype LJ = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        if (i == particleIndex) continue;
        ntype dr = potential.minimalImageDistance(position, particles[i].getPosition());
        // Half the energy to avoid double-counting pairs, correct for potentials
        LJ += potential.lennardJones(dr) * 0.5;
    }
    particles[particleIndex].setPotential(LJ);
    particles[particleIndex].computeKinetic();
    particles[particleIndex].updateEnergy();
}

void particleEnsemble::ensembleStep() {
    for (int i = 0; i < numParticles; ++i) {
        particleEnergyStep(i);
    }
}

void particleEnsemble::storeEnsemble() {
    for (int i = 0; i < numParticles; ++i) {
        particles[i].store();
    }
}

void particleEnsemble::ensembleAverageEnergy() {
    ntype avg = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        avg += particles[i].getEnergy();
    }
    ensembleEnergy = avg / numParticles;
}

void particleEnsemble::Euler(ntype dt, int numSteps, std::ofstream &outFile) {
    ProgressBar pbar(numSteps);

    for (int i = 0; i < numSteps; ++i) {
        ensembleSnapshot(outFile, false);
        storeEnsemble();
        for (int j = 0; j < numParticles; ++j) {
            const Vector oldPos = particles[j].getOldPosition();
            const Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0); 
            for (int k = 0; k < numParticles; ++k) {
                if (k != j) {
                    // Remove the *0.5 in force calculations
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition());
                }
            }

            Vector newAcc = Force / particles[j].getMass();
            Vector newPos = oldPos + oldVel * dt;
            Vector newVel = oldVel + newAcc * dt;

            applyPeriodicBoundary(newPos);
            particles[j].setPosition(newPos);
            particles[j].setVelocity(newVel);
        }
        ensembleStep();
        pbar.update(i);
    }
    pbar.finish();
    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
}

void particleEnsemble::EulerCromer(ntype dt, int numSteps, std::ofstream &outFile) {
    ProgressBar pbar(numSteps);

    for (int i = 0; i < numSteps; ++i) {
        ensembleSnapshot(outFile, false);
        storeEnsemble();
        for (int j = 0; j < numParticles; ++j) {
            const Vector oldPos = particles[j].getOldPosition();
            const Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0); 
            for (int k = 0; k < numParticles; ++k) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition());
                }
            }
            Vector newAcc = Force / particles[j].getMass();

            Vector newVel = oldVel + newAcc * dt;
            Vector newPos = oldPos + newVel * dt;

            applyPeriodicBoundary(newPos);
            particles[j].setPosition(newPos);
            particles[j].setVelocity(newVel);
        }
        ensembleStep();
        pbar.update(i);
    }
    pbar.finish();
    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
}

void particleEnsemble::SpeedVerlet(ntype dt, int numSteps, std::ofstream &outFile) {
    ProgressBar pbar(numSteps);

    for (int i = 0; i < numSteps; i++) {
        ensembleSnapshot(outFile, false);
        storeEnsemble();
        // First half step
        for (int j = 0; j < numParticles; j++) {
            const Vector oldPos = particles[j].getOldPosition();
            const Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition());
                }
            }
            Vector newAcc = Force / particles[j].getMass();
            Vector newPos = oldPos + oldVel * dt + newAcc * (0.5 * dt * dt);
            applyPeriodicBoundary(newPos);
            particles[j].setAcceleration(newAcc);
            particles[j].setPosition(newPos);
        }
        // Second half step
        for (int j = 0; j < numParticles; j++) {
            const Vector pos = particles[j].getPosition();
            const Vector oldVel = particles[j].getOldVelocity();
            const Vector oldAcc = particles[j].getAcceleration();

            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(pos, particles[k].getPosition());
                }
            }
            Vector newAcc = Force / particles[j].getMass();
            Vector newVel = oldVel + (oldAcc + newAcc) * (0.5 * dt);
            particles[j].setVelocity(newVel);
        }

        ensembleStep();
        pbar.update(i);
    }
    pbar.finish();
    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
}

// Predictor-Corrector placeholder (no implementation details provided)
void particleEnsemble::PredictorCorrector(ntype dt, int numSteps, std::ofstream &outFile) {
    ProgressBar pbar(numSteps);
    for (int i = 0; i < numSteps; i++) {
        ensembleSnapshot(outFile, false);
        storeEnsemble();
        // Implementation details would go here.
        ensembleStep();
        pbar.update(i);
    }
}

void particleEnsemble::ThermoEulerCromer(ntype dt, int numSteps, std::ofstream &outFile, ntype temperature, int thermoSteps) {
    ProgressBar pbar(numSteps);
    ntype potentialEnergy = 0.0;
    std::cout << "Executing ThermoEulerCromer\n";   
    for (int i = 0; i < numSteps; ++i) {
        ensembleSnapshot(outFile, false);
        if (i % thermoSteps == 0 && i > 0) {
            // Rescale velocities to maintain temperature
            ntype currentTemperature = 0.0;
            for (int j = 0; j < numParticles; j++) {
                currentTemperature += particles[j].getKinetic() * 2.0 / (3.0 * particles[j].getMass());
            }
            currentTemperature /= numParticles;
            if (currentTemperature > 0.0 && !std::isnan(currentTemperature)) {
                ntype scaling = std::sqrt(temperature / currentTemperature);
                if (!std::isnan(scaling) && !std::isinf(scaling)) {
                    for (int j = 0; j < numParticles; j++) {
                        particles[j].setVelocity(particles[j].getVelocity() * scaling);
                    }
                } else {
                    std::cerr << "[DEBUG] Invalid scaling factor. Skipping velocity scaling.\n";
                }
            } else {
                std::cerr << "[DEBUG] Invalid current temperature. Skipping velocity scaling.\n";
            }
        }
        storeEnsemble();
        for (int j = 0; j < numParticles; ++j) {
            const Vector oldPos = particles[j].getOldPosition();
            const Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0); 
            for (int k = 0; k < numParticles; ++k) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition());
                }
            }
            Vector newAcc = Force / particles[j].getMass();

            Vector newVel = oldVel + newAcc * dt;
            Vector newPos = oldPos + newVel * dt;

            applyPeriodicBoundary(newPos);
            particles[j].setPosition(newPos);
            particles[j].setVelocity(newVel);
        }
        ntype avgPotential = 0.0;
        for (int ip = 0; ip < numParticles; ++ip) {
            avgPotential += particles[ip].getPotential();
        }
        potentialEnergy += avgPotential / numParticles;
        ensembleStep();
        pbar.update(i);
    }
    pbar.finish();
    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
    std::cout << "Average potential energy: " << potentialEnergy / numSteps << std::endl;
}

void particleEnsemble::ThermoSpeedVerlet(ntype dt, int numSteps, std::ofstream &outFile, ntype targetTemperature, int thermoSteps) {
    ProgressBar pbar(numSteps);

    ntype potentialEnergy = 0.0;
    int stabilizationSteps = 200;

    for (int step = 0; step < numSteps; step++) {
        ensembleSnapshot(outFile, false);

        // Velocity rescaling thermostat
        if (step % thermoSteps == 0 && step > 0) {
            ntype currentTemperature = 0.0;
            for (int j = 0; j < numParticles; j++) {
                currentTemperature += particles[j].getKinetic() * 2.0 / (3.0 * particles[j].getMass());
            }
            currentTemperature /= numParticles;
            if (!std::isnan(currentTemperature) && currentTemperature > 0.0) {
                ntype scaling = std::sqrt(targetTemperature / currentTemperature);
                if (!std::isnan(scaling) && !std::isinf(scaling)) {
                    for (int j = 0; j < numParticles; j++) {
                        particles[j].setVelocity(particles[j].getVelocity() * scaling);
                    }
                } else {
                    std::cerr << "[DEBUG] Warning: scaling factor is invalid. Skipping velocity scaling.\n";
                }
            } else {
                std::cerr << "[DEBUG] Warning: currentTemperature is invalid. Skipping velocity scaling.\n";
            }
        }

        storeEnsemble();
        // First half-step (position update)
        for (int j = 0; j < numParticles; j++) {
            const Vector oldPos = particles[j].getOldPosition();
            const Vector oldVel = particles[j].getOldVelocity();

            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition());
                }
            }

            Vector newAcc = Force / particles[j].getMass();

            if (std::isnan(newAcc(0)) || std::isnan(newAcc(1)) || std::isnan(newAcc(2))) {
                std::cerr << "[DEBUG] Warning: Acceleration for particle " << j 
                          << " is NaN.\n";
            }

            Vector newPos = oldPos + oldVel * dt + newAcc * (0.5 * dt * dt);
            applyPeriodicBoundary(newPos);
            particles[j].setAcceleration(newAcc);
            particles[j].setPosition(newPos);
        }

        // Second half-step (velocity update)
        for (int j = 0; j < numParticles; j++) {
            const Vector pos = particles[j].getPosition();
            const Vector oldVel = particles[j].getOldVelocity();
            const Vector oldAcc = particles[j].getAcceleration();

            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(pos, particles[k].getPosition());
                }
            }

            Vector newAcc = Force / particles[j].getMass();
            Vector newVel = oldVel + (oldAcc + newAcc) * (0.5 * dt);
            particles[j].setVelocity(newVel);
        }

        // Update energies
        ensembleStep();

        ntype avgPotential = 0.0;
        for (int j = 0; j < numParticles; j++) {
            avgPotential += particles[j].getPotential();
        }
        avgPotential /= numParticles;

        if (step > stabilizationSteps) {
            potentialEnergy += avgPotential;
        }

        pbar.update(step);
    }
    pbar.finish();

    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
    std::cout << "Average potential energy: " 
              << potentialEnergy / (numSteps - stabilizationSteps) << std::endl;
}
