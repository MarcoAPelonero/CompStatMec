#ifndef INTERACTIONPOTENTIAL_HPP
#define INTERACTIONPOTENTIAL_HPP

#include "vec.hpp"
#include <iostream>
#include <cstdlib>

class interactionPotential {
    protected:
        ntype epsilon;
        ntype sigma;
        ntype L;
        ntype rcut;

    public:
        interactionPotential();
        interactionPotential(ntype boxLenght);
        interactionPotential(ntype sig, ntype eps, ntype boxLenght);
        interactionPotential(ntype sig, ntype eps, ntype boxLenght, ntype rc);
        ~interactionPotential();

        ntype getEpsilon(); 
        ntype getSigma();
        ntype getL();
        ntype getRcut();

        Vector minimalImageDisplacement(const Vector &r1, const Vector &r2);
        ntype minimalImageDistance(const Vector &r1, const Vector &r2);

        ntype lennardJones(ntype r);
        ntype cutLennardJones(ntype r);
        ntype coulomb(ntype r);

        ntype computeForceMagnitudeLennardJones(ntype r);
        ntype computeForceLennardJones(Vector r1, Vector r2);
        ntype computeForceMagnitudeCoulomb(ntype r);
        ntype computeForceCoulomb(Vector r1, Vector r2);
};

#endif // INTERACTIONPOTENTIALS_HPP