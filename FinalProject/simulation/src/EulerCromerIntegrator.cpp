#include "EulerCromerIntegrator.hpp"

void EulerCromerIntegrator::step(std::vector<Particle> &particles, 
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

        velocity += acceleration * dt;
        position += velocity * dt;

        particle.setVelocity(velocity);
        particle.setPosition(position);

        boundaryConditions.apply(particle);
    }
}
