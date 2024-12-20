#ifndef VELOCITY_VERLET_INTEGRATOR_HPP
#define VELOCITY_VERLET_INTEGRATOR_HPP

#include "Integrator.hpp"

// Velocity Verlet Integrator
class SpeedVerletIntegrator : public Integrator {
public:
    void step(std::vector<Particle> &particles, 
              ForceCalculator &forceCalculator,
              BoundaryConditions &boundaryConditions, 
              double dt, 
              double boxLength) override;
};

#endif // VELOCITY_VERLET_INTEGRATOR_HPP
