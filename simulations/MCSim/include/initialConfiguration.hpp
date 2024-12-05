#ifndef INITIALCONFIGURATION_HPP
#define INITIALCONFIGURATION_HPP

#include "generals.hpp"
#include "configuration.hpp"
#include "particles.hpp"
#include <cmath>

class initialConfiguration {
    private: 
        particleEnsemble ensemble;
    public:
        initialConfiguration();
        initialConfiguration(int N, ntype L);
        ~initialConfiguration();

        void initializeSquareLattice(ntype L);
        //  void initializeRandom(ntype L);
        // void initializeRandom(ntype L, ntype T);
        // void initializeRandom(ntype L, ntype T, ntype v);
        // void initializeRandom(ntype L, ntype T, ntype v, ntype a);

        particleEnsemble returnParticleEnsemble();
};

#endif // INITIALCONFIGURATION_HPP