#ifndef PARTICLE_ENSEMBLE_HPP
#define PARTICLE_ENSEMBLE_HPP

#include "particle.hpp"
#include <vector>
#include <fstream>
#include <cmath>
#include <random>

class ParticleEnsemble {
private:
    std::vector<Particle> particles;
    int numParticles;
    double boxLength;
    std::mt19937 rng;
    std::normal_distribution<double> gaussian;
    double totalEnergy, totalKineticEnergy, totalVirialEnergy, totalPotentialEnergy;

    // RDF-related variables
    int rdfBins;
    double rdfCutoff;
    double rdfBinWidth;
    std::vector<double> rdfHistogram;

public:
    ParticleEnsemble(int N, double L)
    : numParticles(N),
      boxLength(L),
      rng(42),
      gaussian(0.0, 1.0),
      totalEnergy(0),
      totalKineticEnergy(0),
      totalVirialEnergy(0),
      totalPotentialEnergy(0),
      rdfBins(1000),       
      rdfCutoff(L/2.0),    
      rdfBinWidth((L/2.0)/100) 
    {
        rdfHistogram.resize(rdfBins, 0.0);
        initializeParticles();
    }

    void initializeParticles();

    std::vector<Particle>& getParticles();
    const std::vector<Particle>& getParticles() const;
    double getBoxLength() const { return boxLength; }
    double getTotalEnergy() const { return totalEnergy; }
    double getTotalKineticEnergy() const { return totalKineticEnergy; }
    double getTotalPotentialEnergy() const { return totalPotentialEnergy; }

    void setParticles(const std::vector<Particle> &p) { particles = p; }
    void setBoxLength(double l) { boxLength = l; }
    void setTotalEnergy(double e) { totalEnergy = e; }
    void setTotalKineticEnergy(double e) { totalKineticEnergy = e; }
    void setTotalPotentialEnergy(double e) { totalPotentialEnergy = e; }

    void addParticle(const Particle &p) { particles.push_back(p); }

    void showParticles() const;

    void periodicBoundaryConditions(Particle &p) const;
    Vec minimumImageConvention(const Particle &p1, const Particle &p2) const;

    double LennardJonesPotential(const Particle &p1, const Particle &p2) const;
    Vec LennardJonesForce(const Particle &p1, const Particle &p2) const;

    void computePotentials();
    void computeKinetics();
    void computeEnergies();
    void computeEnsembleEnergies();
    void computeForces();
    void computeVirial();

    void ensembleThermalize(double temperature, double dt, double tauT = 0.1);
    void ensemblePressurize(double pressure, double dt, double tauP = 0.1);

    void ensembleStepEuler(double dt);
    void ensembleStepEulerCromer(double dt);
    void ensembleStepSpeedVerlet(double dt);

    void ensembleSnapshot(std::ofstream &file);

    // RDF methods
    void resetRDFHistogram();
    void computeRadialDistributionFunctionDirect();
    void computeRadialDistributionFunctionCellMethod();
    void printRadialDistributionFunction(std::ofstream &file);

private:
    struct Cell {
        std::vector<int> particleIndices; 
    };

    void buildCells(std::vector<Cell> &cells, int &numCellsPerDim, double &cellSize) const;
};

#endif // PARTICLE_ENSEMBLE_HPP
