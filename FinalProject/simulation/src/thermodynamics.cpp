#include "thermodynamics.hpp"
#include <cmath>

// Compute Kinetic Energy for All Particles
void Thermodynamics::computeKinetics(std::vector<Particle> &particles) {
    for(auto &p : particles) {
        p.computeKineticEnergy();
    }
}

// Compute Virial from Forces
void Thermodynamics::computeVirial(std::vector<Particle> &particles, double boxLenght, ForceCalculator &forceCalculator) {
    // Reset virial energies
    for(auto &p : particles) {
        p.setVirialEnergy(0.0);
    }

    // Assuming ForceCalculator has a method to compute virial
    forceCalculator.computeVirial(particles, boxLenght);
}

// Compute Total Energies
void Thermodynamics::computeEnergies(double &totalEnergy, double &totalKineticEnergy, 
                                     double &totalVirialEnergy, double &totalPotentialEnergy, 
                                     std::vector<Particle> &particles)
{
    totalVirialEnergy = 0.0;
    totalKineticEnergy = 0.0;
    totalPotentialEnergy = 0.0;
    totalEnergy = 0.0;

    for(auto &p : particles) {
        p.computeTotalEnergy();
        totalVirialEnergy += p.getVirialEnergy();
        totalKineticEnergy += p.getKineticEnergy();
        totalPotentialEnergy += p.getPotentialEnergy();
        totalEnergy += p.getTotalEnergy();
    }
}

// Thermalize Ensemble to Target Temperature
void Thermodynamics::thermalize(std::vector<Particle> &particles, double targetTemperature, double dt, double tauT) {
    double currentTemperature = 0.0;
    for(const auto &p : particles) {
        currentTemperature += p.getKineticEnergy();
    }
    currentTemperature = (2.0 / 3.0) * currentTemperature / particles.size();

    double lambda = std::sqrt(1.0 + (dt / tauT) * (targetTemperature / currentTemperature - 1.0));

    for(auto &p : particles) {
        Vec velocity = p.getVelocity();
        velocity *= lambda;
        p.setVelocity(velocity);
    }

    // Recompute Kinetic Energies
    computeKinetics(particles);
}

// Pressurize Ensemble to Target Pressure
void Thermodynamics::pressurize(std::vector<Particle> &particles, double targetPressure, double dt, double tauP, 
                                 double &boxLength, BoundaryConditions &boundaryConditions) 
{
    double volume = std::pow(boxLength, 3);
    double totalKineticEnergy = 0.0;
    double totalVirialEnergy = 0.0;
    for(const auto &p : particles) {
        totalKineticEnergy += p.getKineticEnergy();
        totalVirialEnergy += p.getVirialEnergy();
    }
    double currentPressure = (2.0 * totalKineticEnergy + totalVirialEnergy) / (3.0 * volume);

    double lambda = 1.0 - (dt / tauP) * (targetPressure - currentPressure);

    // Clamp lambda to prevent excessive scaling
    const double lambda_min = 0.9;
    const double lambda_max = 1.1;
    lambda = std::max(lambda_min, std::min(lambda_max, lambda));

    double boxScale = std::pow(lambda, 1.0 / 3.0);
    boxLength *= boxScale;

    for(auto &p : particles) {
        Vec position = p.getPosition();
        position *= boxScale;
        p.setPosition(position);
        boundaryConditions.apply(p);
    }
}
