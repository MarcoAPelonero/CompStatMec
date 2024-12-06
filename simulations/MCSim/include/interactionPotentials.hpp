#ifndef INTERACTIONPOTENTIALS_HPP
#define INTERACTIONPOTENTIALS_HPP

#include "vec.hpp"

class interactionPotential {
    protected:
        ntype epsilon;
        ntype sigma;
        ntype L;

    public:
        interactionPotential();
        interactionPotential(ntype boxLenght);
        interactionPotential(ntype sig, ntype eps, ntype boxLenght);
        // interactionPotential(ntype eps, ntype sig, ntype rut);
        ~interactionPotential();

        ntype computeDistance(Vector r1, Vector r2);
        ntype lennardJones(ntype r);
};

#endif // INTERACTIONPOTENTIALS_HPP