#ifndef PARTICLEENSEMBLE_HPP
#define PARTICLEENSEMBLE_HPP

#include "particle.hpp"
#include <vector>

class particleEnsemble {
    private: 
        std::vector<Particle> particles;
        int numParticles;
        ntype boxSize;

    public: 
        particleEnsemble(int numParticles, ntype boxSize);
        ~particleEnsemble();


        void show();


};

#endif // PARTICLEENSEMBLE_HPP