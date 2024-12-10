#include "particle.hpp"

Particle::Particle() {
    r = Vector();
    v = Vector();
    a = Vector();
    m = 1.0;
}

Particle::Particle(Vector r) {
    this->r = r;
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

Particle::Particle(const Particle& other) { // Copy constructor
    this->r = other.r;
    this->rold = other.rold;
    this->v = other.v;
    this->vold = other.vold;
    this->a = other.a;
    this->aold = other.aold;
    this->m = other.m;
}

Particle::~Particle() {}

Vector Particle::getPosition() {
    return r; 
}

Vector Particle::getOldPosition() {
    return rold;
}

Vector Particle::getVelocity() {
    return v;
}   

Vector Particle::getOldVelocity() {
    return vold;
}

Vector Particle::getAcceleration() {
    return a;
}

Vector Particle::getOldAcceleration() {
    return aold;
}

ntype Particle::getMass() {
    return m;
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

void Particle::setOldPosition(Vector oldPos) {
    rold = oldPos;
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
