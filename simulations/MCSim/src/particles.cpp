#include "particles.hpp"

Particle::Particle() {
    r = Vector();
    v = Vector();
    a = Vector();
    m = 1.0;
}

Particle::Particle(Vector r) {
    this->v = r;
    v = Vector();
    a = Vector();
    m = 1.0;
}

Particle::Particle(Vector r, Vector v, Vector a, ntype m) {
    this->r = r;
    this->v = v;
    this->a = a;
    this->m = m;
}

Particle::~Particle() {}


Vector Particle::getPosition() {
    return r; 
}

void Particle::store() {
    rold = r;
}

void Particle::restore() {
    r = rold;
}

void Particle::random(ntype L) {
    r.random(L);
}

particleEnsemble::particleEnsemble(int N) {
    particles = std::vector<Particle>(N);
}

 particleEnsemble::~particleEnsemble() {}

void particleEnsemble::initializeRandom(ntype L) {
    for (auto& p : particles) {
        p.random(L);
    }
}


Particle& particleEnsemble::operator()(int i) {
    return particles[i];
}

void particleEnsemble::store() {
    for (auto& p : particles) {
        p.store();
    }
}

void particleEnsemble::restore() {
    for (auto& p : particles) {
        p.restore();
    }
}