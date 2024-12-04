#include "monteCarloMove.hpp"

MonteCarloMove::MonteCarloMove() {
    ensemble = particleEnsemble();
    T = 1.0;
    v = 1.0;
    a = 1.0;

    p = Particle();
    numParticles = particleEnsemble().getNumParticles();
}

MonteCarloMove::MonteCarloMove(particleEnsemble ensemble, ntype T) {
    this->ensemble = ensemble;
    this->T = T;
    v = 1.0;
    a = 1.0;

    numParticles = particleEnsemble().getNumParticles();
    particleExtractor.setDistribution(0, ensemble.getNumParticles() - 1);
}

MonteCarloMove::MonteCarloMove(particleEnsemble ensemble, ntype T, ntype v, ntype a) {
    this->ensemble = ensemble;
    this->T = T;
    this->v = v;
    this->a = a;

    particleExtractor.setDistribution(0, ensemble.getNumParticles() - 1);
}

MonteCarloMove::~MonteCarloMove() {}


void MonteCarloMove::chooseParticle() {
    particleIndex = particleExtractor();
    p = ensemble(particleIndex);
}

void MonteCarloMove::savePosition() {
    p.store();
}

void MonteCarloMove::moveParticle() {
    // Example implementation: Randomly move each particle
    Vector randOrient, displacement; 
    randOrient.random();
    displacement = randOrient * delta;
    p.setPosition(p.getPosition() + displacement);
}

void MonteCarloMove::computeEnergy() {
    ensemble.calculateEnergy();
}

void MonteCarloMove::computeEnergyDifference() {
    ntype deltaE = 0.0;
    Vector r_i_new = p.getPosition();
    Vector r_i_old = p.getOldPosition();

    for (int j = 0; j < ensemble.getNumParticles(); ++j) {
        if (j != particleIndex) {
            Particle& pj = ensemble(j);
            Vector r_j = pj.getPosition();

            ntype dr_old = computeDistance(r_i_new, r_j);
            ntype dr_new = computeDistance(r_i_old, r_j);

            ntype V_new = lennardJones(dr_new);
            ntype V_old = lennardJones(dr_old);

            deltaE += V_new - V_old;
        }
    }

    // Store deltaE in a class member variable for later use
    this->deltaE = deltaE;
}

void MonteCarloMove::acceptMove() {
    ntype acceptanceProbability = exp(-deltaE / T);

    if (deltaE < 0 || rng.ranf() < acceptanceProbability) {
        ensemble(particleIndex) = p;
        ensemble.updateEnergy(deltaE);
    } else {
        p.restore();
    }
}

void MonteCarloMove::MetropolisStep() {
    chooseParticle();
    savePosition();
    moveParticle();
    computeEnergyDifference();
    acceptMove();
}