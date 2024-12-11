#include "particleEnsemble.hpp"
#include <cmath>
#include <random>
#include <iostream>

particleEnsemble::particleEnsemble(int numParticles, ntype boxSize)
    : numParticles(numParticles), boxSize(boxSize) {
    int particlesPerSide = std::ceil(std::cbrt(numParticles));
    ntype spacing = boxSize / particlesPerSide;

    std::default_random_engine generator;
    std::normal_distribution<ntype> distribution(0.0, 0.1); // smaller std to reduce initial explosion

    int count = 0;
    for (int i = 0; i < particlesPerSide && count < numParticles; ++i) {
        for (int j = 0; j < particlesPerSide && count < numParticles; ++j) {
            for (int k = 0; k < particlesPerSide && count < numParticles; ++k) {
                Particle p(Vector{i * spacing, j * spacing, k * spacing});

                // Smaller initial velocities:
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

void particleEnsemble::show() {
    for (int i = 0; i < numParticles; ++i) {
        std::cout << "Particle " << i+1 << ":\n";
        particles[i].show();
    }
}

void particleEnsemble::stepEuler(int i, ntype dt) {
    Vector r1 = particles[i].getPosition();
    Vector v  = particles[i].getVelocity();

    Vector fij{0.0, 0.0, 0.0};
    ntype potentialEnergy = 0.0;

    for (int j = 0; j < numParticles; ++j) {
        if (i == j) continue;

        Vector dr = potential.minimalImageDisplacement(r1, particles[j].getOldPosition());
        ntype r = dr.modulus();
        ntype pe = potential.lennardJones(r);
        ntype f = potential.computeForceMagnitude(r);
        potentialEnergy += pe;
        fij += (dr / r) * f;
    }

    Vector a = fij / particles[i].getMass();
    Vector newPosition = r1 + v * dt;
    Vector newVelocity = v + a * dt;

    // Periodic boundary conditions
    for (int d = 0; d < dim; ++d) {
        if (newPosition(d) < 0) {
            newPosition(d) += boxSize;
        } else if (newPosition(d) >= boxSize) {
            newPosition(d) -= boxSize;
        }
    }

    particles[i].setPosition(newPosition);
    particles[i].setVelocity(newVelocity);
    particles[i].updateEnergy(potentialEnergy);
}

void particleEnsemble::stepEulerCromer(int i, ntype dt) {
    Vector r1 = particles[i].getPosition();
    Vector v  = particles[i].getVelocity();

    Vector fij{0.0, 0.0, 0.0};
    ntype potentialEnergy = 0.0;

    for (int j = 0; j < numParticles; ++j) {
        if (i == j) continue;

        Vector dr = potential.minimalImageDisplacement(r1, particles[j].getOldPosition());
        ntype r = dr.modulus();
        ntype pe = potential.lennardJones(r);
        ntype f = potential.computeForceMagnitude(r);
        potentialEnergy += pe;
        fij += (dr / r) * f;
    }

    Vector a = fij / particles[i].getMass();
    Vector newVelocity = v + a * dt;
    Vector newPosition = r1 + newVelocity * dt;

    // Periodic boundary conditions
    for (int d = 0; d < dim; ++d) {
        if (newPosition(d) < 0) {
            newPosition(d) += boxSize;
        } else if (newPosition(d) >= boxSize) {
            newPosition(d) -= boxSize;
        }
    }

    particles[i].setPosition(newPosition);
    particles[i].setVelocity(newVelocity);
    particles[i].updateEnergy(potentialEnergy);
}

void particleEnsemble::updateEnsemble() {
    for (int i = 0; i < numParticles; ++i) {
        particles[i].store();
    }
}

void particleEnsemble::ensembleStep(ntype dt) {
    for (int i = 0; i < numParticles; ++i) {
        stepEuler(i, dt);
    }
    updateEnsemble();
}

void particleEnsemble::ensembleSnapshot(std::ofstream &outFile, bool writeHeader) {
    if (writeHeader) {
        outFile << "x y z vx vy vz energy\n";
    }
    for (int i = 0; i < numParticles; ++i) {
        Vector pos = particles[i].getPosition();
        Vector vel = particles[i].getVelocity();
        ntype E = particles[i].getEnergy();
        outFile << pos(0) << " " << pos(1) << " " << pos(2) << " "
                << vel(0) << " " << vel(1) << " " << vel(2) << " "
                << E << "\n";
    }
    outFile << "\n";
}
