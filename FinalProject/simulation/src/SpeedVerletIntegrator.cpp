#include "SpeedVerletIntegrator.hpp"

void SpeedVerletIntegrator::step(std::vector<Particle> &particles, 
                                     ForceCalculator &forceCalculator,
                                     BoundaryConditions &boundaryConditions, 
                                     double dt, 
                                     double boxLength)
{
    forceCalculator.computeForces(particles, boxLength);

    std::vector<Vec> oldAccelerations;
    oldAccelerations.reserve(particles.size());

    for(auto &particle : particles) {
        Vec acceleration = particle.getForce() / particle.getMass();
        oldAccelerations.push_back(acceleration);

        Vec position = particle.getPosition();
        Vec velocity = particle.getVelocity();
        position += velocity * dt + acceleration * 0.5 * dt * dt;
        particle.setPosition(position);

        boundaryConditions.apply(particle);
    }

    forceCalculator.computeForces(particles, boxLength);

    for(size_t i = 0; i < particles.size(); ++i) {
        Vec newAcceleration = particles[i].getForce() / particles[i].getMass();
        Vec velocity = particles[i].getVelocity();
        velocity += 0.5 * (oldAccelerations[i] + newAcceleration) * dt;
        particles[i].setVelocity(velocity);
    }
}
