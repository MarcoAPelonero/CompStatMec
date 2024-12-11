#include "monteCarloMove.hpp"

MonteCarloMove::MonteCarloMove(particleEnsemble& ensemble, ntype T, ntype delta)
    : ensemble(ensemble), T(T), 
    delta(delta), L(ensemble.getBoxLength()), 
    deltaE(0.0), numParticles(ensemble.getNumParticles()), 
    particleIndex(0), acceptedMoves(0), rejectedMoves(0) {}

MonteCarloMove::~MonteCarloMove() {}

MonteCarloMove& MonteCarloMove::operator=(const MonteCarloMove& other) {
    if (this == &other) return *this;
    ensemble = other.ensemble;
    T = other.T;
    delta = other.delta;
    L = other.L;
    deltaE = other.deltaE;
    numParticles = other.numParticles;
    particleIndex = other.particleIndex;
    acceptedMoves = other.acceptedMoves;
    rejectedMoves = other.rejectedMoves;
    return *this;
}

void MonteCarloMove::MetropolisStep() {
    selectParticle();
    Vector dr = calculateDisplacement();
    moveParticle(dr);
    deltaE = ensemble.calculateEnergyDifference(particleIndex);

    // std::cout << "Delta E: " << deltaE << std::endl;
    // std::cout << "New Enegy: " << newEnergy << std::endl;
    // std::cout << "Old Energy: " << oldEnergy << std::endl;
    
    checkAcceptance();
}

void MonteCarloMove::selectParticle() {
    particleIndex = particleExtractor() % numParticles;
    // std::cout << "Selected particle: " << particleIndex << std::endl;
}

Vector MonteCarloMove::calculateDisplacement() {
    Vector dr;
    for (int i = 0; i < dim; ++i) {
        dr(i) = rng.ranf() * delta - delta / 2.0;
    }
    return dr;
}

void MonteCarloMove::moveParticle(const Vector& dr) {
    Particle& p = ensemble(particleIndex);
    ensemble.store();
    Vector r = p.getPosition();
    Vector newR = r + dr;
    p.setPosition(newR);
}

void MonteCarloMove::checkAcceptance() {
    ntype r = rng.ranf();
    ntype boltzmann = std::exp(-deltaE / T);
    if (r < boltzmann) {
        acceptedMoves++;
        ensemble.updateEnsemble(particleIndex, deltaE);
    } else {
        ensemble.restoreEnergy();
        ensemble.restoreParticle(particleIndex);
        rejectedMoves++;
    }
}

void MonteCarloMove::printAcceptanceRatio() const {
    int totalMoves = acceptedMoves + rejectedMoves;
    double acceptanceRatio = static_cast<double>(acceptedMoves) / totalMoves;
    std::cout << "Acceptance ratio: " << acceptanceRatio << std::endl;
}