#include "particles.hpp"

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


particleEnsemble::particleEnsemble() {
    particles = std::vector<Particle>(0);
    numParticles = 0;
}
particleEnsemble::particleEnsemble(int N, ntype L) {
    particles = std::vector<Particle>(N);
    numParticles = N;
    boxLength = L;
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
    
    for (int i = 0; i < numParticles; ++i) {
        for (int j = 0; j < i; ++j) { // Changed j loop to run from 0 to i-1
            Vector r_i = particles[i].getPosition();
            Vector r_j = particles[j].getPosition();
            ntype distance = computeDistance(r_i, r_j);
            energy += lennardJones(distance);
        }
    }
    return energy;
}

ntype particleEnsemble::getEnergy() {
    return ensembleEnergy;
}

void particleEnsemble::setInitialEnergy() {
    ensembleEnergy = calculateEnergy();
}

// void particleEnsemble::setEnergy(ntype newEnergy) {
//     ensembleEnergy = newEnergy;
// }

void particleEnsemble::updateEnergy(ntype deltaEnergy) {
    ensembleEnergy += deltaEnergy; 
}

void particleEnsemble::saveSnapshotPosition(std::ofstream& file, int step) {
    if (!file.is_open()) {
        std::cerr << "Error: File is not open for writing." << std::endl;
        return;
    }

    // Write the step number
    file << "Step " << step << std::endl;

    // Write the positions of all particles
    for (auto& p : particles) {
        Vector pos = p.getPosition();
        for (int i = 0; i < dim; ++i) {
            file << pos.get(i) << " ";
        }
        file << std::endl;
    }
}

void particleEnsemble::saveSnapshotVelocity(std::ofstream& file) {
    if (!file.is_open()) {
        std::cerr << "Error: File is not open for writing." << std::endl;
        return;
    }

    for (auto& p : particles) {
        Vector vel = p.getVelocity();
        for (int i = 0; i < dim; ++i) {
            file << vel(i) << " ";
        }
        file << "\n";
    }
}

void particleEnsemble::ensembleSnapshot(std::ofstream& file) {
    if (!file.is_open()) {
        std::cerr << "Error: File is not open for writing." << std::endl;
        return;
    }

    for (auto& p : particles) {
        Vector pos = p.getPosition();
        Vector vel = p.getVelocity();
        for (int i = 0; i < dim; ++i) {
            file << pos.get(i) << " ";
        }
        for (int i = 0; i < dim; ++i) {
            file << vel.get(i) << " ";
        }
        file << "\n";
    }
}