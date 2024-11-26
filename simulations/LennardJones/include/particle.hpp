#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Vector.hpp"
#include "config.hpp"

class Particle {
protected:
    Vector position;
    Vector old_position;
    ntype sigma, epsilon, rcut;

public:
    Particle() : sigma(0), epsilon(0), rcut(0) {}
    Particle(const Vector& pos, ntype s, ntype e, ntype rc)
        : position(pos), sigma(s), epsilon(e), rcut(rc) {}

    void setPosition(const Vector& pos) {
        position = pos;
    }

    void setSigma(ntype s) {
        sigma = s;
    }

    void setEpsilon(ntype e) {
        epsilon = e;
    }

    void setRcut(ntype rc) {
        rcut = rc;
    }

    Vector getPosition() const {
        return position;
    }

    ntype getSigma() const {
        return sigma;
    }

    ntype getEpsilon() const {
        return epsilon;
    }

    ntype getRcut() const {
        return rcut;
    }

    void store() {
        old_position = position;
    }

    void restore() {
        position = old_position;
    }
};

class LennardJonesParticle : public Particle {

public:
    LennardJonesParticle(const Vector& pos, ntype s, ntype e, ntype rc)
        : Particle(pos, s, e, rc) {}

    ntype vij(const LennardJonesParticle& P, Vector& L)
    {
        ntype ene;
        Vector Dr;

        Dr = position - P.position;

        // MINIMUM IMAGE CONVENTION 
        Dr = Dr - L.mulcw(rint(divcw(Dr,L))); // Dr - L * rint(Dr / L)

        ntype rn = Dr.norm();
        ntype rsq = rn * rn;

        if (rsq < rcut * rcut) // Use 'rcut' from base class
            ene = 4.0 * epsilon * (pow(sigma / rn, 12.0) - pow(sigma / rn, 6.0));
        else
            ene = 0.0;

        return ene;
    }
};

class FENEParticle : public Particle {
private:
    ntype k;  // Spring constant
    ntype R0; // Maximum extension

public:
    FENEParticle() : k(0), R0(0) {}
    FENEParticle(const Vector& pos, ntype s, ntype e, ntype rc, ntype k_, ntype R0_)
        : Particle(pos, s, e, rc), k(k_), R0(R0_) {}

    void setK(ntype k_) {
        k = k_;
    }

    void setR0(ntype R0_) {
        R0 = R0_;
    }

    ntype getK() const {
        return k;
    }

    ntype getR0() const {
        return R0;
    }
};

#endif // PARTICLE_HPP