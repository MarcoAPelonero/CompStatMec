#ifndef MONTECARLOMOVE_HPP
#define MONTECARLOMOVE_HPP

#include "particles.hpp"
#include "configuration.hpp"
#include "randNumGen.hpp"
#include "interactionPotentials.hpp"
#include "vec.hpp"

class MonteCarloMove: public interactionPotential {
    private:
        particleEnsemble ensemble;
        ntype T;
        ntype v;
        ntype a;

        Particle p;
        ntype deltaE;
        int numParticles;
        int particleIndex;
    public:
        MonteCarloMove();
        MonteCarloMove(particleEnsemble ensemble, ntype T);
        MonteCarloMove(particleEnsemble ensemble, ntype T, ntype v, ntype a);
        ~MonteCarloMove();

        void MetropolisStep();

        void chooseParticle();
        void moveParticle();
        void savePosition();
        void computeEnergy();
        void computeEnergyDifference();
        void acceptMove();
};

#endif // MONTECARLOMOVE_HPP