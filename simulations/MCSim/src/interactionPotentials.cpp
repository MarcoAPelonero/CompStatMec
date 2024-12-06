#include "interactionPotentials.hpp"


interactionPotential::interactionPotential() {
    sigma = 1.0;
    L = 10.0;
    epsilon = 1.0;
}
// src/interactionPotentials.cpp
interactionPotential::interactionPotential(ntype boxLength)
    : L(boxLength) {
        sigma = 1.0;
        epsilon = 1.0;
    }

interactionPotential::interactionPotential(ntype sig, ntype eps, ntype boxLength) {
    sigma = sig;
    epsilon = eps;
    L = boxLength;
}

//interactionPotential::interactionPotential(ntype eps, ntype sig, ntype rut) {
//    epsilon = eps;
//    sigma = sig;
//    rcut = rc;
//}

interactionPotential::~interactionPotential() {}

ntype interactionPotential::computeDistance(Vector r1, Vector r2) {
    Vector dr = r1 - r2;
    // Apply minimum image convention
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
    // td::cout << "Lennard-Jones potential: " << V << std::endl;
    return V;
}
