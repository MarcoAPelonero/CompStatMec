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

        Vector minimalImageDisplacement(Vector &r1, const Vector &r2);
        ntype minimalImageDistance(Vector &r1, const Vector &r2);

        ntype lennardJones(ntype r);
        ntype cutLennardJones(ntype r);

        ntype computeForceMagnitude(ntype r);
        ntype computeForce(Vector r1, Vector r2);
        ntype computeForceMagnittude(ntype r);  
};

#endif // INTERACTIONPOTENTIALS_HPP