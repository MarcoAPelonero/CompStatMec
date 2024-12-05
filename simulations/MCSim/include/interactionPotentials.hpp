#ifndef INTERACTIONPOTENTIALS_HPP
#define INTERACTIONPOTENTIALS_HPP

#include "vec.hpp"

class interactionPotential {
    protected:
        ntype epsilon;
        ntype sigma;

    public:
        interactionPotential();
        interactionPotential(ntype sig);
        // interactionPotential(ntype eps, ntype sig, ntype rut);
        ~interactionPotential();

        ntype computeDistance(Vector r1, Vector r2);
        ntype lennardJones(ntype r);
};

#endif // INTERACTIONPOTENTIALS_HPP