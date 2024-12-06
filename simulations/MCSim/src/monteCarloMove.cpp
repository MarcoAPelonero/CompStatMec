#include "monteCarloMove.hpp"

MonteCarloMove::MonteCarloMove(particleEnsemble& ensemble, ntype T, ntype L)
    : interactionPotential(L), ensemble(ensemble), T(T), L(L), 
      acceptedMoves(0), rejectedMoves(0) {
    particleExtractor.setDistribution(0, ensemble.getNumParticles() - 1);
    particleIndex = 0;
    numParticles = ensemble.getNumParticles();
    this->L = L;
}

MonteCarloMove::MonteCarloMove(particleEnsemble& ensemble, ntype T, ntype v, ntype a) 
    : interactionPotential(), ensemble(ensemble), T(T), v(v), a(a), 
      acceptedMoves(0), rejectedMoves(0) {
    numParticles = ensemble.getNumParticles();
    particleExtractor.setDistribution(0, ensemble.getNumParticles() - 1);
    particleIndex = 0;
}


MonteCarloMove::~MonteCarloMove() {}

void MonteCarloMove::chooseParticle() {
    particleIndex = particleExtractor();
}

void MonteCarloMove::moveParticle() {
    Vector randOrient, displacement, pos, newPos;
    pos = ensemble(particleIndex).getPosition();
    // std::cout << std::endl;
    // pos.show("Current position: ");
    randOrient = Vector::randomVector();
    displacement = randOrient;
    newPos = pos + displacement;
    // newPos.show("New position: ");
    for (int i = 0; i < dim; ++i) {
        if (newPos(i) >= L) newPos(i) -= L;
        if (newPos(i) < 0) newPos(i) += L;
    }

    ensemble(particleIndex).setPosition(newPos);
}

void MonteCarloMove::computeEnergyDifference() {
    ntype deltaE = 0.0;
    Vector r_i_new = ensemble(particleIndex).getPosition();
    Vector r_i_old = ensemble(particleIndex).getOldPosition();

    for (int j = 0; j < numParticles; ++j) {
        if (j != particleIndex) {
            Particle& pj = ensemble(j);
            Vector r_j = pj.getPosition();

            ntype dr_old = computeDistance(r_i_old, r_j);
            ntype dr_new = computeDistance(r_i_new, r_j);

            ntype V_new = lennardJones(dr_new);
            ntype V_old = lennardJones(dr_old);

            deltaE += V_new - V_old;
        }
    }

    this->deltaE = deltaE;
}

void MonteCarloMove::acceptMove() {
    ntype acceptanceProbability = exp(-deltaE / T);
    ntype r = rng.ranf();
    // std::cout << "Accepting move: " << acceptanceProbability << std::endl;
    // std::cout << "Random number: " << r << std::endl;
    if (r < acceptanceProbability) {
        // Move accepted
        acceptedMoves++;
        ensemble(particleIndex).setOldPosition(ensemble(particleIndex).getPosition());
        ensemble.updateEnergy(deltaE);
    } else {
        // Move rejected
        rejectedMoves++;
        ensemble(particleIndex).restore();
    }
}

void MonteCarloMove::MetropolisStep() {
    chooseParticle();
    ensemble(particleIndex).store(); // Store the current position
    moveParticle();
    computeEnergyDifference();
    acceptMove();
}

void MonteCarloMove::printAcceptanceRatio() {
    std::cout << std::endl << "Acceptance ratio: "
              << static_cast<double>(acceptedMoves) / (acceptedMoves + rejectedMoves)
              << std::endl;
}