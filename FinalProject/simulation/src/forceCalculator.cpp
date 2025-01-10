#include "ForceCalculator.hpp"
#include <cmath>

// Calculate displacement using Minimum Image Convention
Vec ForceCalculator::minimumImage(const Particle &p1, const Particle &p2, double boxLength) const {
    Vec r = p1.getPosition() - p2.getPosition();
    double halfBox = boxLength / 2.0;

    // Apply Minimum Image Convention
    if (r.getX() > halfBox) r.setX(r.getX() - boxLength);
    if (r.getX() < -halfBox) r.setX(r.getX() + boxLength);
    if (r.getY() > halfBox) r.setY(r.getY() - boxLength);
    if (r.getY() < -halfBox) r.setY(r.getY() + boxLength);
    if (r.getZ() > halfBox) r.setZ(r.getZ() - boxLength);
    if (r.getZ() < -halfBox) r.setZ(r.getZ() + boxLength);

    return r;
}

// Compute Lennard-Jones Potential between two particles
double ForceCalculator::LennardJonesPotential(const Particle &p1, const Particle &p2, double boxLength) const {
    Vec r = minimumImage(p1, p2, boxLength);
    double r_len = r.length();

    if (r_len == 0.0 || r_len > 4.0) return 0.0;

    double r2 = r_len * r_len;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    return 4.0 * (1.0 / r12 - 1.0 / r6);
}

// Compute Lennard-Jones Force between two particles
Vec ForceCalculator::LennardJonesForce(const Particle &p1, const Particle &p2, double boxLength) const {
    Vec r = minimumImage(p1, p2, boxLength);
    double r_len = r.length();

    if (r_len == 0.0 || r_len > 4.0) return Vec(0.0, 0.0, 0.0);

    double r2 = r_len * r_len;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    double force_magnitude = 24.0 * (2.0 / r12 - 1.0 / r6) / r2;

    return r * force_magnitude;
}

// Compute Potentials for All Particle Pairs
void ForceCalculator::computePotentials(std::vector<Particle> &particles, double boxLength) {
    // Reset potential energies
    for(auto &particle : particles) {
        particle.setPotentialEnergy(0.0);
    }

    // Compute pairwise potentials
    for(size_t i = 0; i < particles.size(); ++i) {
        for(size_t j = i + 1; j < particles.size(); ++j) {
            double potential = LennardJonesPotential(particles[i], particles[j], boxLength);
            particles[i].setPotentialEnergy(particles[i].getPotentialEnergy() + 0.5 * potential);
            particles[j].setPotentialEnergy(particles[j].getPotentialEnergy() + 0.5 * potential);
        }
    }
}

// Compute Forces for All Particle Pairs
void ForceCalculator::computeForces(std::vector<Particle> &particles, double boxLength) {
    // Reset forces
    for(auto &particle : particles) {
        particle.setForce(Vec(0.0, 0.0, 0.0));
    }

    // Compute pairwise forces
    for(size_t i = 0; i < particles.size(); ++i) {
        for(size_t j = i + 1; j < particles.size(); ++j) {
            Vec force = LennardJonesForce(particles[i], particles[j], boxLength);
            particles[i].addForce(force);
            particles[j].addForce(force * -1.0);
        }
    }
}

// Compute Virial Contributions and Update Particles
void ForceCalculator::computeVirial(std::vector<Particle> &particles, double boxLength) {
    // Reset virial energies
    for(auto &particle : particles) {
        particle.setVirialEnergy(0.0);
    }

    // Compute pairwise virial contributions
    for(size_t i = 0; i < particles.size(); ++i) {
        for(size_t j = i + 1; j < particles.size(); ++j) {
            Vec r = minimumImage(particles[i], particles[j], boxLength);
            Vec force = LennardJonesForce(particles[i], particles[j], boxLength);
            double virial = r.dot(force);

            // Each pair contributes equally to both particles
            particles[i].setVirialEnergy(particles[i].getVirialEnergy() + 0.5 * virial);
            particles[j].setVirialEnergy(particles[j].getVirialEnergy() + 0.5 * virial);
        }
    }
}
