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
    double totalEnergy, totalKineticEnergy, totalPotentialEnergy;

public:
    ParticleEnsemble(int N, double L)
        : numParticles (N), boxLength(L), rng(42), gaussian(0.0, 1.0), totalEnergy(0), totalKineticEnergy(0), totalPotentialEnergy(0) {
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

    void ensembleThermalize(double temperature);
    void ensembleStepEuler(double dt);  
    void ensembleStepEulerCromer(double dt);

    void ensembleSnapshot(std::ofstream &file);
};

#endif // PARTICLE_ENSEMBLE_HPP
