#ifndef EULER_CROMER_INTEGRATOR_HPP
#define EULER_CROMER_INTEGRATOR_HPP

#include "Integrator.hpp"

// Euler-Cromer Integrator
class EulerCromerIntegrator : public Integrator {
public:
    void step(std::vector<Particle> &particles, 
              ForceCalculator &forceCalculator,
              BoundaryConditions &boundaryConditions, 
              double dt, 
              double boxLength) override;
};

#endif // EULER_CROMER_INTEGRATOR_HPP
