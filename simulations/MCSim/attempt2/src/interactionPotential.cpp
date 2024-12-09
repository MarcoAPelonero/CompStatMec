#include "interactionPotential.hpp"


interactionPotential::interactionPotential() {
    sigma = 1.0;
    L = 10.0;
    epsilon = 1.0;
    rcut = 2.5;
}

interactionPotential::interactionPotential(ntype boxLength)
    : L(boxLength) {
        sigma = 1.0;
        epsilon = 1.0;
        rcut = 2.5;
    }

interactionPotential::interactionPotential(ntype sig, ntype eps, ntype boxLength) {
    sigma = sig;
    epsilon = eps;
    L = boxLength;
    rcut = 2.5;
}

interactionPotential::interactionPotential(ntype sig, ntype eps, ntype boxLength, ntype rc) {
    sigma = sig;
    epsilon = eps;
    L = boxLength;
    rcut = rc;
}

interactionPotential::~interactionPotential() {}

ntype interactionPotential::getEpsilon() {
    return epsilon;
}

ntype interactionPotential::getSigma() {
    return sigma;
}

ntype interactionPotential::getL() {
    return L;
}

ntype interactionPotential::getRcut() {
    return rcut;
}

ntype interactionPotential::computeDistance(Vector r1, Vector r2) {
    Vector dr = r1 - r2;
    for (int i = 0; i < dim; ++i) {
        if (dr(i) > L / 2.0) {
            dr(i) -= L;
        } else if (dr(i) < -L / 2.0) {
            dr(i) += L;
        }
    }  
    double r = dr.modulus();
    return r;
}

ntype interactionPotential::lennardJones(ntype r) {  
    ntype V = 4.0 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
    return V;
}

ntype interactionPotential::cutLennardJones(ntype r) {
    if (r < rcut) {
        return 4.0 * (std::pow(1.0 / r, 12) - std::pow(1.0 / r, 6));
    } else {
        return 0.0;
    }
}