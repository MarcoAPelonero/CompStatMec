#ifndef EULER_INTEGRATOR_HPP
#define EULER_INTEGRATOR_HPP

#include "Integrator.hpp"

// Euler Integrator
class EulerIntegrator : public Integrator {
public:
    void step(std::vector<Particle> &particles, 
              ForceCalculator &forceCalculator,
              BoundaryConditions &boundaryConditions, 
              double dt, 
              double boxLength) override;
};

#endif // EULER_INTEGRATOR_HPP
