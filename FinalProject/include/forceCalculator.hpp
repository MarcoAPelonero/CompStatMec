#ifndef FORCE_CALCULATOR_HPP
#define FORCE_CALCULATOR_HPP

#include "particle.hpp"
#include "Vec.hpp"

#include <vector>

class ForceCalculator {
public:
    ForceCalculator() = default;
    ~ForceCalculator() = default;

    // Compute potentials for all particles
    void computePotentials(std::vector<Particle> &particles, double boxLength);

    // Compute forces and apply to particles
    void computeForces(std::vector<Particle> &particles, double boxLength);

    // Compute virial contributions and update particles
    void computeVirial(std::vector<Particle> &particles, double boxLength);

private:
    // Lennard-Jones Potential
    double LennardJonesPotential(const Particle &p1, const Particle &p2, double boxLength) const;

    // Lennard-Jones Force
    Vec LennardJonesForce(const Particle &p1, const Particle &p2, double boxLength) const;

    // Calculate the displacement vector using Minimum Image Convention
    Vec minimumImage(const Particle &p1, const Particle &p2, double boxLength) const;
};

#endif // FORCE_CALCULATOR_HPP
