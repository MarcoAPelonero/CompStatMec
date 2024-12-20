#ifndef THERMODYNAMICS_HPP
#define THERMODYNAMICS_HPP

#include "particle.hpp"
#include "forceCalculator.hpp"
#include "boundaryConditions.hpp"

#include <vector>

class Thermodynamics {
public:
    Thermodynamics() = default;
    ~Thermodynamics() = default;

    // Compute Kinetic Energy
    void computeKinetics(std::vector<Particle> &particles);

    // Compute Virial
    void computeVirial(std::vector<Particle> &particles, double boxLenght, ForceCalculator &forceCalculator);

    // Compute Total Energies
    void computeEnergies(double &totalEnergy, double &totalKineticEnergy, 
                         double &totalVirialEnergy, double &totalPotentialEnergy, 
                         std::vector<Particle> &particles);

    // Thermalization and Pressurization
    void thermalize(std::vector<Particle> &particles, double targetTemperature, double dt, double tauT);
    void pressurize(std::vector<Particle> &particles, double targetPressure, double dt, double tauP, 
                   double &boxLength, BoundaryConditions &boundaryConditions);
};

#endif // THERMODYNAMICS_HPP
