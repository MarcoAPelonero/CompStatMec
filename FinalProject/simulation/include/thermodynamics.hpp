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

    void computeKinetics(std::vector<Particle> &particles);

    void computeVirial(std::vector<Particle> &particles, double boxLenght, ForceCalculator &forceCalculator);

    void computeEnergies(double &totalEnergy, double &totalKineticEnergy, 
                         double &totalVirialEnergy, double &totalPotentialEnergy, 
                         std::vector<Particle> &particles);

    void thermalize(std::vector<Particle> &particles, double targetTemperature, double dt, double tauT);
    void pressurize(std::vector<Particle> &particles, double targetPressure, double dt, double tauP, 
                   double &boxLength, BoundaryConditions &boundaryConditions);
};

#endif // THERMODYNAMICS_HPP
