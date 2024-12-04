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

Vector Particle::getOldPosition() {
    return rold;
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

void Particle::setPosition(Vector pos) {
    r = pos;
}

void Particle::setVelocity(Vector vel) {
    v = vel;
}

Particle& Particle::operator=(Particle p) {
    if (this == &p)
        return *this;

    this->r = p.r;
    this->rold = p.rold;
    this->v = p.v;
    this->vold = p.vold;
    this->a = p.a;
    this->aold = p.aold;
    this->m = p.m;

    return *this;
}


particleEnsemble::particleEnsemble() {
    particles = std::vector<Particle>(0);
    numParticles = 0;
}
particleEnsemble::particleEnsemble(int N) {
    particles = std::vector<Particle>(N);
    numParticles = N;
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

int particleEnsemble::getNumParticles() {
    return numParticles;
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

void particleEnsemble::initParticle(int i, Vector v) {
    particles[i].setPosition(v);
}

ntype particleEnsemble::calculateEnergy() {
    ntype energy = 0;
    // Placeholder for potential calculation
    // Replace this with the actual potential calculation once defined
    for (int i = 0; i < numParticles; ++i) {
        for (int j = i + 1; j < numParticles; ++j) {
            // energy += potential(particles[i], particles[j]);
        }
    }
    ensembleEnergy = energy;
    return energy;
}

ntype particleEnsemble::getEnergy() {
    return ensembleEnergy;
}

void particleEnsemble::updateEnergy(ntype newEnergy) {
    ensembleEnergy = newEnergy;
}