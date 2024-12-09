#ifndef PARTICLEENSEMBLE_HPP
#define PARTICLEENSEMBLE_HPP

#include "vec.hpp"
#include "particle.hpp"
#include "interactionPotential.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>

class particleEnsemble {
private:
    std::vector<Particle> particles;
    int numParticles;  
    ntype ensembleEnergy;
    ntype pastEnergy;

    ntype boxLength;

    interactionPotential potential; // Interaction potential member
public:
    // Constructors
    particleEnsemble();
    particleEnsemble(int N, ntype L, double sigma =1.0, double epsilon =1.0);
    ~particleEnsemble();

    void initializeEnsemble();
    void initializeRandom(ntype L);

    Particle& operator()(int i);
    void initParticle(int i, Vector v);

    void store();
    void restore();

    int getNumParticles() const;
    ntype getEnergy() const;
    void setEnergy(ntype newEnergy);
    void setOldEnergy(ntype oldEnergy);
    void restoreEnergy();
    
    ntype calculateEnergy();
    ntype calculateEnergyDifference(int particleIndex);  
    void updateEnergy(ntype deltaEnergy);
    void updateParticle(int i, Vector v);
    void saveSnapshotPosition(std::ofstream& file, int step);

    ntype getBoxLength() const;

    interactionPotential& getPotential() { return potential; }
    const interactionPotential& getPotential() const { return potential; }

    void clearParticles();
    void addParticle(const Particle& p);

    Vector applyMinimumImageConvention(const Vector& position) const;
    void moveParticle(int i, const Vector& dr);

private:
    int integerCubeRoot(int N) const;
};

#endif // PARTICLEENSEMBLE_HPP
