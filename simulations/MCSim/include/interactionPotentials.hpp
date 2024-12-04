#ifndef INTERACTIONPOTENTIALS_HPP
#define INTERACTIONPOTENTIALS_HPP

#include "vec.hpp"
#include "configuration.hpp"

class interactionPotential {
    private:
        ntype epsilon;
        ntype sigma;
        ntype rcut;
    public:
        interactionPotential();
        interactionPotential(ntype sig);
        interactionPotential(ntype eps, ntype sig, ntype rc);
        ~interactionPotential();

        ntype computeDistance(Vector r1, Vector r2);
        ntype lennardJones(ntype r);
};

#endif // INTERACTIONPOTENTIALS_HPP