#include "interactionPotentials.hpp"

interactionPotential::interactionPotential() {
    epsilon = 1.0;
    sigma = 1;
}

interactionPotential::interactionPotential(ntype sig) {
    sigma = sig;
}

//interactionPotential::interactionPotential(ntype eps, ntype sig, ntype rut) {
//    epsilon = eps;
//    sigma = sig;
//    rcut = rc;
//}

interactionPotential::~interactionPotential() {}

ntype interactionPotential::computeDistance(Vector r1, Vector r2) {
    Vector dr = r1 - r2;
    // std::cout << "Position 1: ";
    //r1.show();
    // std::cout << "Position 2: ";
    // r2.show();
    // std::cout << "Difference vector: ";
    // dr.show();
    ntype r = dr.modulus();
    // std::cout << "Distance r: " << r << std::endl;
    return r;
}

ntype interactionPotential::lennardJones(ntype r) {  
    if (r>rc) {
        return 0;
    }
    ntype V = 4.0 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
    // td::cout << "Lennard-Jones potential: " << V << std::endl;
    return V;
}
