#ifndef PARTICLEENSEMBLE_HPP
#define PARTICLEENSEMBLE_HPP

#include "particle.hpp"
#include "interactionPotential.hpp"
#include <vector>
#include <fstream>

class particleEnsemble {
    private:
        std::vector<Particle> particles;
        int numParticles;
        ntype boxSize;
        interactionPotential potential;

    public:
        particleEnsemble(int numParticles, ntype boxSize);
        ~particleEnsemble();

        Particle& operator()(int i);
        int getNumParticles();

        void show();
        void stepEuler(int i, ntype dt);
        void stepEulerCromer(int i, ntype dt);
        void updateEnsemble();
        void ensembleStep(ntype dt);

        void ensembleSnapshot(std::ofstream &outFile, bool writeHeader = false);
};

#endif // PARTICLEENSEMBLE_HPP
