#include "particleEnsemble.hpp"
#include <cmath>
#include <random>
#include <iostream>

particleEnsemble::particleEnsemble(int numParticles, ntype boxSize)
    : numParticles(numParticles), boxSize(boxSize), potential(boxSize) {
    int particlesPerSide = std::ceil(std::cbrt(numParticles));
    ntype spacing = boxSize / particlesPerSide;

    std::default_random_engine generator;
    std::normal_distribution<ntype> distribution(0.0, 0.5); // smaller std to reduce initial explosion

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

void particleEnsemble::updateEnsemble() {
    for (int i = 0; i < numParticles; ++i) {
        particles[i].store();
    }
}

// Helper method to compute force and potential energy for a single particle i at a given position
void particleEnsemble::computeForceAndEnergyForParticle(int i, Vector &pos, Vector &force, ntype &potentialEnergy) {
    force = Vector{0.0, 0.0, 0.0};
    potentialEnergy = 0.0;

    for (int j = 0; j < numParticles; ++j) {
        if (i == j) continue;
        Vector dr = potential.minimalImageDisplacement(pos, particles[j].getOldPosition());
        ntype r = dr.modulus();
        ntype pe = potential.lennardJones(r);
        ntype f = potential.computeForceMagnitude(r);
        potentialEnergy += pe;
        force += (dr / r) * f;
    }
}

// Apply periodic boundary conditions
inline void particleEnsemble::applyPeriodicBoundary(Vector &pos) {
    for (int d = 0; d < dim; ++d) {
        if (pos(d) < 0) {
            pos(d) += boxSize;
        } else if (pos(d) >= boxSize) {
            pos(d) -= boxSize;
        }
    }
}

// Euler step for a single particle
void particleEnsemble::stepEuler(int i, ntype dt) {
    Vector r1 = particles[i].getPosition();
    Vector v  = particles[i].getVelocity();

    Vector fij;
    ntype potentialEnergy;
    computeForceAndEnergyForParticle(i, r1, fij, potentialEnergy);

    Vector a = fij / particles[i].getMass();
    Vector newPosition = r1 + v * dt;
    Vector newVelocity = v + a * dt;

    applyPeriodicBoundary(newPosition);

    particles[i].setPosition(newPosition);
    particles[i].setVelocity(newVelocity);
    particles[i].updateEnergy(potentialEnergy);
}

// Euler-Cromer step for a single particle
void particleEnsemble::stepEulerCromer(int i, ntype dt) {
    Vector r1 = particles[i].getOldPosition();
    Vector v  = particles[i].getOldVelocity();

    Vector fij;
    ntype potentialEnergy;
    computeForceAndEnergyForParticle(i, r1, fij, potentialEnergy);

    Vector a = fij / particles[i].getMass();
    Vector newVelocity = v + a * dt;
    Vector newPosition = r1 + newVelocity * dt;

    applyPeriodicBoundary(newPosition);

    particles[i].setPosition(newPosition);
    particles[i].setVelocity(newVelocity);
    particles[i].updateEnergy(potentialEnergy);
}

//NOT USABLE AT THE MOMENT
void particleEnsemble::stepLeapFrog(int i, ntype dt) {
    // Positions and velocities from previous time
    Vector r_old = particles[i].getOldPosition();
    Vector v_half_old = particles[i].getOldVelocity(); // assumed half-step velocity

    // Update position
    Vector r_new = r_old + v_half_old * dt;
    applyPeriodicBoundary(r_new);

    // Compute new acceleration at t+dt
    Vector fij_new;
    ntype potentialEnergy;
    computeForceAndEnergyForParticle(i, r_new, fij_new, potentialEnergy);
    Vector a_new = fij_new / particles[i].getMass();

    // Update half-step velocity: v(t+3dt/2) = v(t+dt/2) + a(t+dt)*dt
    Vector v_half_new = v_half_old + a_new * dt;

    particles[i].setPosition(r_new);
    particles[i].setVelocity(v_half_new);
    particles[i].updateEnergy(potentialEnergy);
}

// Velocity Verlet step for a single particle
void particleEnsemble::stepVelocityVerlet(int i, ntype dt) {
    Vector r1 = particles[i].getOldPosition();
    Vector v  = particles[i].getOldVelocity();

    Vector fij;
    ntype potentialEnergy;
    computeForceAndEnergyForParticle(i, r1, fij, potentialEnergy);
    Vector a = fij / particles[i].getMass();

    // Position half-step
    Vector newPosition = r1 + (v * dt) + a * (0.5 * dt * dt);
    applyPeriodicBoundary(newPosition);

    // Compute force with updated position
    Vector fij_new;
    ntype dummyEnergy; // potentialEnergy stored for the initial pos only, or reassign here if desired
    computeForceAndEnergyForParticle(i, newPosition, fij_new, dummyEnergy);
    Vector a_new = fij_new / particles[i].getMass();

    // Velocity full-step
    Vector newVelocity = v + (a + a_new) * (0.5 * dt);

    particles[i].setPosition(newPosition);
    particles[i].setVelocity(newVelocity);
    particles[i].updateEnergy(potentialEnergy);
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
