#include "EulerIntegrator.hpp"

void EulerIntegrator::step(std::vector<Particle> &particles, 
                           ForceCalculator &forceCalculator,
                           BoundaryConditions &boundaryConditions, 
                           double dt, 
                           double boxLength)
{
    // Compute Forces
    forceCalculator.computeForces(particles, boxLength);

    // Update Positions and Velocities
    for(auto &particle : particles) {
        Vec acceleration = particle.getForce() / particle.getMass();
        Vec velocity = particle.getVelocity();
        Vec position = particle.getPosition();

        position += velocity * dt;
        velocity += acceleration * dt;

        particle.setPosition(position);
        particle.setVelocity(velocity);

        // Apply Boundary Conditions
        boundaryConditions.apply(particle);
    }
}
