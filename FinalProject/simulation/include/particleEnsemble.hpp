#ifndef PARTICLE_ENSEMBLE_HPP
#define PARTICLE_ENSEMBLE_HPP

#include "EulerIntegrator.hpp"
#include "EulerCromerIntegrator.hpp"
#include "SpeedVerletIntegrator.hpp"
#include <iostream>
#include <cmath>

#include "particle.hpp"
#include "ForceCalculator.hpp"
#include "Integrator.hpp"
#include "RDFCalculator.hpp"
#include "BoundaryConditions.hpp"
#include "Thermodynamics.hpp"
#include "model.hpp"
#include "json.hpp"

#include <vector>
#include <fstream>
#include <random>

class ParticleEnsemble {
private:
    std::vector<Particle> particles;
    int numParticles;
    double boxLength;
    std::mt19937 rng;
    std::normal_distribution<double> gaussian;

    double totalEnergy, totalKineticEnergy, totalVirialEnergy, totalPotentialEnergy;

    ForceCalculator forceCalculator;
    Thermodynamics thermodynamics;
    RDFCalculator rdfCalculator;
    BoundaryConditions boundaryConditions;

    RDFDensityModel mlModel;
public:
    ParticleEnsemble(int N, double L);

    void initializeParticles();

    std::vector<Particle>& getParticles();
    const std::vector<Particle>& getParticles() const;
    double getBoxLength() const;
    double getTotalEnergy() const;
    double getTotalKineticEnergy() const;
    double getTotalPotentialEnergy() const;

    void setParticles(const std::vector<Particle> &p);
    void setBoxLength(double l);
    void setTotalEnergy(double e);
    void setTotalKineticEnergy(double e);
    void setTotalPotentialEnergy(double e);
    void addParticle(const Particle &p);

    void showParticles() const;

    void ensembleThermalize(double temperature, double dt, double tauT = 0.1);
    void ensemblePressurize(double pressure, double dt, double tauP = 0.1);

    void ensembleStepEuler(double dt);
    void ensembleStepEulerCromer(double dt);
    void ensembleStepSpeedVerlet(double dt);

    void ensembleSnapshot(std::ofstream &file);

    void resetRDFHistogram();
    void computeRadialDistributionFunctionDirect();
    void computeRadialDistributionFunctionOptimized();
    void computeThermodynamicsFromRDF(std::ofstream &file);
    double computeEnergyPerParticleFromRDF();
    double computeThermodynamicsFromML();
    std::pair<double, double> computeThermodynamicsFromPythonML();
    void printRadialDistributionFunction(std::ofstream &file);
};

#endif // PARTICLE_ENSEMBLE_HPP