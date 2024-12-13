#include "interactionPotential.hpp"
#include <cmath>

// Constructors remain the same

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

ntype interactionPotential::getEpsilon() { return epsilon; }
ntype interactionPotential::getSigma() { return sigma; }
ntype interactionPotential::getL() { return L; }
ntype interactionPotential::getRcut() { return rcut; }

// Minimal image displacement using L/2.0
Vector interactionPotential::minimalImageDisplacement(const Vector &r1, const Vector &r2) {
    Vector dr = r1 - r2;
    for (int i = 0; i < dim; ++i) {
        if (dr(i) > L/2) {
            dr(i) -= L;
        } else if (dr(i) < -L/2) {
            dr(i) += L;
        }
    }
    return dr;
}

ntype interactionPotential::minimalImageDistance(const Vector &r1, const Vector &r2) {
    Vector dr = minimalImageDisplacement(r1, r2);
    return dr.modulus();
}

ntype interactionPotential::lennardJones(ntype r) {
    ntype sr6 = std::pow(sigma / r, 6);
    ntype sr12 = sr6 * sr6;
    return 4.0 * epsilon * (sr12 - sr6);
}

ntype interactionPotential::cutLennardJones(ntype r) {
    if (r < rcut) {
        ntype sr6 = std::pow(sigma / r, 6);
        ntype sr12 = sr6 * sr6;
        ntype lj = 4.0 * epsilon * (sr12 - sr6);
        ntype sr6c = std::pow(sigma / rcut, 6);
        ntype sr12c = sr6c * sr6c;
        ntype lj_c = 4.0 * epsilon * (sr12c - sr6c);
        return lj - lj_c;
    } else {
        return 0.0;
    }
}

// Compute only the magnitude of the force, given r
ntype interactionPotential::computeForceMagnitudeLennardJones(ntype r) {
    if (r <= 0) return 0.0;
    ntype sr = sigma / r;
    ntype sr6 = std::pow(sr, 6);
    ntype sr12 = sr6 * sr6;
    // f(r) = 24ε (2sr^12 - sr^6)/r
    ntype f = 24.0 * epsilon * (2.0 * sr12 - sr6) / r;
    return f;
}

// Deprecated: If needed, ensure this also uses minimal image properly.
Vector interactionPotential::computeForceLennardJones(Vector r1, Vector r2) {
    Vector dr = minimalImageDisplacement(r1, r2);
    ntype r = dr.modulus();
    return dr * computeForceMagnitudeLennardJones(r);
}

ntype interactionPotential::coulomb(ntype r) {
    if (r <= 0) return 0.0;
    return 1.0 / r;
}

ntype interactionPotential::computeForceMagnitudeCoulomb(ntype r) {
    if (r <= 0) return 0.0;
    // f(r) = 1 / r^2
    return 1.0 / (r * r);
}

// Deprecated: If needed, ensure this also uses minimal image properly.
ntype interactionPotential::computeForceCoulomb(Vector r1, Vector r2) {
    Vector dr = minimalImageDisplacement(r1, r2);
    ntype r = dr.modulus();
    return computeForceMagnitudeCoulomb(r);
}