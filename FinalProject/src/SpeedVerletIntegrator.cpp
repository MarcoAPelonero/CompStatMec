#include "SpeedVerletIntegrator.hpp"

void SpeedVerletIntegrator::step(std::vector<Particle> &particles, 
                                     ForceCalculator &forceCalculator,
                                     BoundaryConditions &boundaryConditions, 
                                     double dt, 
                                     double boxLength)
{
    // Compute Initial Forces
    forceCalculator.computeForces(particles, boxLength);

    // Store Current Accelerations
    std::vector<Vec> oldAccelerations;
    oldAccelerations.reserve(particles.size());

    for(auto &particle : particles) {
        Vec acceleration = particle.getForce() / particle.getMass();
        oldAccelerations.push_back(acceleration);

        // Update Positions
        Vec position = particle.getPosition();
        Vec velocity = particle.getVelocity();
        position += velocity * dt + acceleration * 0.5 * dt * dt;
        particle.setPosition(position);

        // Apply Boundary Conditions
        boundaryConditions.apply(particle);
    }

    // Compute New Forces
    forceCalculator.computeForces(particles, boxLength);

    // Update Velocities with Average of Old and New Accelerations
    for(size_t i = 0; i < particles.size(); ++i) {
        Vec newAcceleration = particles[i].getForce() / particles[i].getMass();
        Vec velocity = particles[i].getVelocity();
        velocity += 0.5 * (oldAccelerations[i] + newAcceleration) * dt;
        particles[i].setVelocity(velocity);
    }
}
