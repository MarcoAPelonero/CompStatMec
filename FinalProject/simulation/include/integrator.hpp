#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "particle.hpp"
#include "forceCalculator.hpp"
#include "boundaryConditions.hpp"

#include <vector>

// Abstract Base Class for Integrators
class Integrator {
public:
    virtual ~Integrator() = default;

    // Pure virtual method to perform a simulation step
    virtual void step(std::vector<Particle> &particles, 
                      ForceCalculator &forceCalculator,
                      BoundaryConditions &boundaryConditions, 
                      double dt, 
                      double boxLength) = 0;
};

#endif // INTEGRATOR_HPP
