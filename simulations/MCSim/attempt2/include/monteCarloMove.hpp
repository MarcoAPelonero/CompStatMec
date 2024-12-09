// MonteCarloMove.hpp
#ifndef MONTECARLOMOVE_HPP
#define MONTECARLOMOVE_HPP

#include "particleEnsemble.hpp"
#include "interactionPotential.hpp"
#include "randNumGen.hpp"
#include <iostream>

class MonteCarloMove {
private:
    particleEnsemble& ensemble; 
    ntype T;             // Temperature
    ntype delta;         // Maximum displacement
    ntype L;             // Box length

    ntype deltaE;        // Energy difference for the current move

    int numParticles;    // Total number of particles
    int particleIndex;   // Index of the selected particle

    int acceptedMoves;   // Counter for accepted moves
    int rejectedMoves;   // Counter for rejected moves

public:
    MonteCarloMove(particleEnsemble& ensemble, ntype T, ntype delta);
    ~MonteCarloMove();

    Vector calculateDisplacement();
    void moveParticle(const Vector& dr);
    void selectParticle();
    void checkAcceptance();

    void MetropolisStep();
    void printAcceptanceRatio() const;

    MonteCarloMove& operator=(const MonteCarloMove& other);

    int getAcceptedMoves() const { return acceptedMoves; }
    int getRejectedMoves() const { return rejectedMoves; }
    double getAcceptanceRatio() const { return static_cast<double>(acceptedMoves) / (acceptedMoves + rejectedMoves); }

    void incrementAcceptedMoves() { acceptedMoves++; } 
    void incrementRejectedMoves() { rejectedMoves++; }

    ntype getTemperature() const { return T; }
    ntype getDelta() const { return delta; }
    int getParticleIndex() const { return particleIndex; }

    void setDeltaE(ntype newDeltaE) { deltaE = newDeltaE; }
    void setTemperature(ntype newT) { T = newT; }
};

#endif // MONTECARLOMOVE_HPP