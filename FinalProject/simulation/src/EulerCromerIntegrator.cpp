#include "EulerCromerIntegrator.hpp"

void EulerCromerIntegrator::step(std::vector<Particle> &particles, 
                                  ForceCalculator &forceCalculator,
                                  BoundaryConditions &boundaryConditions, 
                                  double dt, 
                                  double boxLength)
{
    // Compute Forces
    forceCalculator.computeForces(particles, boxLength);

    // Update Velocities and Positions
    for(auto &particle : particles) {
        Vec acceleration = particle.getForce() / particle.getMass();
        Vec velocity = particle.getVelocity();
        Vec position = particle.getPosition();

        velocity += acceleration * dt;
        position += velocity * dt;

        particle.setVelocity(velocity);
        particle.setPosition(position);

        // Apply Boundary Conditions
        boundaryConditions.apply(particle);
    }
}
