#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "particle.hpp"
#include "forceCalculator.hpp"
#include "boundaryConditions.hpp"

#include <vector>

class Integrator {
public:
    virtual ~Integrator() = default;

    virtual void step(std::vector<Particle> &particles, 
                      ForceCalculator &forceCalculator,
                      BoundaryConditions &boundaryConditions, 
                      double dt, 
                      double boxLength) = 0;
};

#endif // INTEGRATOR_HPP
