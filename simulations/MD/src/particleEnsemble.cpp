#include "particleEnsemble.hpp"
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <limits>


particleEnsemble::particleEnsemble(int numParticles, ntype boxSize)
    : numParticles(numParticles), boxSize(boxSize), potential(boxSize) {
    int particlesPerSide = std::ceil(std::cbrt(numParticles));
    ntype spacing = boxSize / particlesPerSide;

    std::default_random_engine generator;
    std::normal_distribution<ntype> distribution(0.0, 1); 

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
                if (count >= numParticles) break;
            }
        }
    }

    ntype minDistance = std::numeric_limits<ntype>::max();
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            ntype distance = (particles[i].getPosition() - particles[j].getPosition()).modulus();
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

void particleEnsemble::applyPeriodicBoundary(Vector &pos) {
    for (int i = 0; i < dim; ++i) {
        while (pos(i) < 0) pos(i) += boxSize;
        while (pos(i) >= boxSize) pos(i) -= boxSize;
    }
}

void particleEnsemble::ensembleSnapshot(std::ofstream &outFile, bool writeHeader) {
    if (writeHeader) {
        outFile << "x y z vx vy vz energy\n";
    }
    for (int i = 0; i < numParticles; ++i) {
        Vector pos = particles[i].getPosition();
        Vector vel = particles[i].getVelocity();
        ntype E = particles[i].getPotential();
        outFile << pos(0) << " " << pos(1) << " " << pos(2) << " "
                << vel(0) << " " << vel(1) << " " << vel(2) << " "
                << E << "\n";
    }
    outFile << "\n";
}

void particleEnsemble::particleEnergyStep(int particleIndex){
    Vector position = particles[particleIndex].getPosition();
    ntype LJ = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        if (i == particleIndex) continue;
        ntype dr = potential.minimalImageDistance(position, particles[i].getPosition());
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
            Vector oldPos = particles[j].getOldPosition();
            Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0); 
            for (int k = 0; k < numParticles; ++k) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition()) * 0.5;
                }
            }
            Vector newAcc = Force / particles[j].getMass();

            Vector newPos = oldPos + oldVel * dt;
            Vector newVel = oldVel + newAcc * dt;

            applyPeriodicBoundary(newPos);
            for (int k = 0; k < dim; ++k) {
                if (newPos(k) < 0 || newPos(k) >= boxSize) {
                    std::cerr << "Particle " << j << " escaped the box\n";
                    std::cerr << "Position: ";
                    newPos.show();
                }
            }
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
            Vector oldPos = particles[j].getOldPosition();
            Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0); 
            for (int k = 0; k < numParticles; ++k) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition()) * 0.5;
                }
            }
            Vector newAcc = Force / particles[j].getMass();

            Vector newVel = oldVel + newAcc * dt;
            Vector newPos = oldPos + newVel * dt;

            applyPeriodicBoundary(newPos);
            for (int k = 0; k < dim; ++k) {
                if (newPos(k) < 0 || newPos(k) >= boxSize) {
                    std::cerr << "Particle " << j << " escaped the box\n";
                    std::cerr << "Position: ";
                    newPos.show();
                }
            }
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
        for (int j = 0; j < numParticles; j++) {
            Vector oldPos = particles[j].getOldPosition();
            Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition()) * 0.5;
                }
            }
            Vector newAcc = Force / particles[j].getMass();

            Vector newPos = oldPos + oldVel * dt + newAcc * 0.5 * dt * dt;
            applyPeriodicBoundary(newPos);
            
            particles[j].setAcceleration(newAcc);
            particles[j].setPosition(newPos);
        }
        for (int j = 0; j < numParticles; j++) {
            Vector oldPos = particles[j].getOldPosition();
            Vector pos = particles[j].getPosition();
            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++)
                Force += potential.computeForceLennardJones(pos, particles[k].getPosition()) * 0.5;
            Vector newAcc = Force / particles[j].getMass();
            Vector newVel = particles[j].getOldVelocity() + (particles[j].getAcceleration() + newAcc) * 0.5 * dt;
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

void particleEnsemble::PredictorCorrector(ntype dt, int numSteps, std::ofstream &outFile) {
    ProgressBar pbar(numSteps);

    for (int i = 0; i < numSteps; i++) {
        ensembleSnapshot(outFile, false);
        storeEnsemble();
        for (int j = 0; j < numParticles; j++) {
            Vector r = particles[j].getPosition();
            Vector v = particles[j].getVelocity();

            Vector Force = Vector(0.0, 0.0, 0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(r, particles[k].getPosition()) * 0.5;
                }
            }
            Vector rPredict = r + v * dt;
            
        }
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
            ntype currentTemperature = 0.0;
            for (int j = 0; j < numParticles; j++) {
                currentTemperature += particles[j].getKinetic() * 2.0 / (3.0 * particles[j].getMass());
            }
            // std::cout << "Temp per particle: " << currentTemperature << std::endl;
            currentTemperature /= numParticles;
            // std::cout << "Current temperature: " << currentTemperature << std::endl;
            ntype scaling = std::sqrt(temperature / currentTemperature);
            for (int j = 0; j < numParticles; j++) {
                particles[j].setVelocity(particles[j].getVelocity() * scaling);
            }
        }
        storeEnsemble();
        for (int j = 0; j < numParticles; ++j) {
            Vector oldPos = particles[j].getOldPosition();
            Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0); 
            for (int k = 0; k < numParticles; ++k) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition()) * 0.5;
                }
            }
            Vector newAcc = Force / particles[j].getMass();

            Vector newVel = oldVel + newAcc * dt;
            Vector newPos = oldPos + newVel * dt;

            applyPeriodicBoundary(newPos);
            for (int k = 0; k < dim; ++k) {
                if (newPos(k) < 0 || newPos(k) >= boxSize) {
                    std::cerr << "Particle " << j << " escaped the box\n";
                    std::cerr << "Position: ";
                    newPos.show();
                }
            }
            particles[j].setPosition(newPos);
            particles[j].setVelocity(newVel);
        }
        ntype avgPotential = 0.0;
        for (int i = 0; i < numParticles; ++i) {
            avgPotential += particles[i].getPotential();
        }   
        potentialEnergy += avgPotential / numParticles;
        ensembleStep();
        pbar.update(i);
    }
    pbar.finish();
    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
    std::cout << "Average potential energy: " << potentialEnergy / (numSteps) << std::endl;
}

void particleEnsemble::ThermoSpeedVerlet(ntype dt, int numSteps, std::ofstream &outFile, ntype temperature, int thermoSteps) {
    ProgressBar pbar(numSteps);
    ntype potentialEnergy = 0.0;
    int stabilizationSteps = 1000;
    for (int i = 0; i < numSteps; i++) {
        ensembleSnapshot(outFile, false);
        storeEnsemble();
        if (i % thermoSteps == 0 && i > 0) {
            ntype currentTemperature = 0.0;
            for (int j = 0; j < numParticles; j++) {
                currentTemperature += particles[j].getKinetic() * 2.0 / (3.0 * particles[j].getMass());
            }
            // std::cout << "Temp per particle: " << currentTemperature << std::endl;
            currentTemperature /= numParticles;
            // std::cout << "Current temperature: " << currentTemperature << std::endl;
            ntype scaling = std::sqrt(temperature / currentTemperature);
            for (int j = 0; j < numParticles; j++) {
                particles[j].setVelocity(particles[j].getVelocity() * scaling);
            }
        }
        storeEnsemble();
        for (int j = 0; j < numParticles; j++) {
            Vector oldPos = particles[j].getOldPosition();
            Vector oldVel = particles[j].getOldVelocity();
            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < numParticles; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(oldPos, particles[k].getOldPosition());
                }
            }
            Vector newAcc = Force / (particles[j].getMass());

            Vector newPos = oldPos + oldVel * dt + newAcc * 0.5 * dt * dt;
            applyPeriodicBoundary(newPos);
            
            particles[j].setAcceleration(newAcc);
            particles[j].setPosition(newPos);
        }
        for (int j = 0; j < numParticles; j++) {
            Vector oldPos = particles[j].getOldPosition();
            Vector pos = particles[j].getPosition();
            Vector Force(0.0,0.0,0.0);
            for (int k = 0; k < j; k++) {
                if (k != j) {
                    Force += potential.computeForceLennardJones(pos, particles[k].getPosition());
                }
            }
            Vector newAcc = Force / (particles[j].getMass());
            Vector newVel = particles[j].getOldVelocity() + (particles[j].getAcceleration() + newAcc) * 0.5 * dt;
            particles[j].setVelocity(newVel);
        }
        ntype avgPotential = 0.0;
        for (int i = 0; i < numParticles; ++i) {
            avgPotential += particles[i].getPotential();
        }   
        if (i > stabilizationSteps)
            potentialEnergy += avgPotential / numParticles;
        ensembleStep();
        pbar.update(i);
    }
    pbar.finish();
    ensembleAverageEnergy();
    std::cout << std::setprecision(10);
    std::cout << "Average energy: " << ensembleEnergy << std::endl;
    std::cout << "Average potential energy: " << potentialEnergy / (numSteps - stabilizationSteps) << std::endl;

}