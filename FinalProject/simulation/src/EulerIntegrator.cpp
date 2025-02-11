#include "EulerIntegrator.hpp"

void EulerIntegrator::step(std::vector<Particle> &particles, 
                           ForceCalculator &forceCalculator,
                           BoundaryConditions &boundaryConditions, 
                           double dt, 
                           double boxLength)
{
    forceCalculator.computeForces(particles, boxLength);

    for(auto &particle : particles) {
        Vec acceleration = particle.getForce() / particle.getMass();
        Vec velocity = particle.getVelocity();
        Vec position = particle.getPosition();

        position += velocity * dt;
        velocity += acceleration * dt;

        particle.setPosition(position);
        particle.setVelocity(velocity);

        boundaryConditions.apply(particle);
    }
}
