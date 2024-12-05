#ifndef MONTECARLOMOVE_HPP
#define MONTECARLOMOVE_HPP

#include "particles.hpp"
#include "configuration.hpp"
#include "randNumGen.hpp"
#include "interactionPotentials.hpp"
#include "initialConfiguration.hpp"
#include "vec.hpp"

class MonteCarloMove: public interactionPotential {
    private:
        particleEnsemble& ensemble; 
        ntype T;
        ntype v;
        ntype a;

        ntype L;
        ntype deltaE;

        int numParticles;
        int particleIndex;

        int acceptedMoves;
        int rejectedMoves;

    public:
        MonteCarloMove(particleEnsemble& ensemble, ntype T, ntype L);
        MonteCarloMove(particleEnsemble& ensemble, ntype T, ntype v, ntype a);
        ~MonteCarloMove();

        void MetropolisStep();

        void chooseParticle();
        void moveParticle();
        void computeEnergyDifference();
        void acceptMove();

        void printAcceptanceRatio();
};

#endif // MONTECARLOMOVE_HPP