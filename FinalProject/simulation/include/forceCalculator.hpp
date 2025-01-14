#ifndef FORCE_CALCULATOR_HPP
#define FORCE_CALCULATOR_HPP

#include "particle.hpp"
#include "Vec.hpp"

#include <vector>

class ForceCalculator {
public:
    ForceCalculator() = default;
    ~ForceCalculator() = default;

    void computePotentials(std::vector<Particle> &particles, double boxLength);

    void computeForces(std::vector<Particle> &particles, double boxLength);

    void computeVirial(std::vector<Particle> &particles, double boxLength);

private:
    double LennardJonesPotential(const Particle &p1, const Particle &p2, double boxLength) const;

    Vec LennardJonesForce(const Particle &p1, const Particle &p2, double boxLength) const;

    Vec minimumImage(const Particle &p1, const Particle &p2, double boxLength) const;
};

#endif // FORCE_CALCULATOR_HPP
