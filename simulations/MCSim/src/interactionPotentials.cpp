#include "interactionPotentials.hpp"

interactionPotential::interactionPotential() {
    epsilon = 1.0;
    sigma = 1.0;
    rcut = 2.5;
}

interactionPotential::interactionPotential(ntype sig) {
    sigma = sig;
}

interactionPotential::interactionPotential(ntype eps, ntype sig, ntype rc) {
    epsilon = eps;
    sigma = sig;
    rcut = rc;
}

interactionPotential::~interactionPotential() {}

ntype interactionPotential::computeDistance(Vector r1, Vector r2) {
    Vector dr = r1 - r2;
    ntype r = dr.modulus();
    return r;
}

ntype interactionPotential::lennardJones(ntype r) {
    ntype r2 = r * r;
    ntype r6 = r2 * r2 * r2;
    ntype r12 = r6 * r6;
    ntype V = 4.0 * epsilon * (sigma * sigma / r12 - sigma * sigma / r6);
    return V;
}

