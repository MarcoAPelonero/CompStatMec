#include "particle.hpp"


Particle::Particle() {
    this->r = Vector();
    this->v = Vector();
    this->a = Vector();
    this->m = 1.0;

    this->energy = 0;
}

Particle::Particle(Vector r) {
    this->r = r;
    this->v = Vector();
    this->a = Vector();
    this->m = 1.0;

    this->energy = 0;
}

Particle::Particle(Vector r, Vector v, Vector a, ntype m) {
    this->r = r;
    this->v = v;
    this->a = a;
    this->m = m;

    this->energy = 0;
}

Particle::Particle(const Particle& other) { // Copy constructor
    this->r = other.r;
    this->rold = other.rold;
    this->v = other.v;
    this->vold = other.vold;
    this->a = other.a;
    this->aold = other.aold;
    this->m = other.m;
    this->energy = other.energy;
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

ntype Particle::getEnergy() {
    return energy;
}

void Particle::store() {
    rold = r;
    vold = v;
    aold = a;   
}

void Particle::restore() {
    r = rold;
    v = vold; 
    a = aold;
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

void Particle::show() {
    std::cout << "Position: ";
    r.show();
    std::cout << "Velocity: ";
    v.show();
    std::cout << "Acceleration: ";
    a.show();
    std::cout << "Mass: " << m << std::endl;

    std::cout << "Energy: " << energy << std::endl;
    std::cout << std::endl;
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

void Particle::updateEnergy(ntype pe) {
    ntype kineticEnergy = 0.5 * m * v.modulus() * v.modulus();
    energy = kineticEnergy + pe;
}