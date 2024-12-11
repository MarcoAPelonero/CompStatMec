// particleEnsemble.cpp
#include "particleEnsemble.hpp"
#include <cmath>
#include <cassert>

// Helper function to compute integer cube root
int particleEnsemble::integerCubeRoot(int N) const {
    if (N <= 0) return 0;
    int root = std::round(std::cbrt(N));
    // Ensure root^3 <= N < (root+1)^3
    while (root*root*root > N) root--;
    while ((root+1)*(root+1)*(root+1) <= N) root++;
    return root;
}

// Default constructor
particleEnsemble::particleEnsemble() : numBins(0), numParticles(0), ensembleEnergy(0.0), pastEnergy(0.0), boxLength(0.0), potential() {
    particles = std::vector<Particle>();
    radialDistFuncHisto = std::vector<ntype>();
}

particleEnsemble::particleEnsemble(int N, ntype L, double sigma, double epsilon) 
    : numParticles(N), ensembleEnergy(0.0), pastEnergy(0.0), boxLength(L), potential(sigma, epsilon, L) {
    if (N < 0) {
        throw std::invalid_argument("Number of particles cannot be negative.");
    }
    if (L <= 0.0) {
        throw std::invalid_argument("Box length must be positive.");
    }
    particles = std::vector<Particle>(N);
    radialDistFuncHisto = std::vector<ntype>();
}
// Parameterized constructor
particleEnsemble::particleEnsemble(int numBins, int N, ntype L, double sigma, double epsilon) 
    : numParticles(N), ensembleEnergy(0.0), pastEnergy(0.0), boxLength(L), potential(sigma, epsilon, L) {
    if (N < 0) {
        throw std::invalid_argument("Number of particles cannot be negative.");
    }
    if (L <= 0.0) {
        throw std::invalid_argument("Box length must be positive.");
    }
    particles = std::vector<Particle>(N);
    radialDistFuncHisto = std::vector<ntype>(numBins);
}

// Destructor
particleEnsemble::~particleEnsemble() {}

void particleEnsemble::initializeEnsemble() {
    if (numParticles == 0) return;

    int latticePerSide = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(numParticles))));
    if (latticePerSide < 1) latticePerSide = 1;

    double spacing = boxLength / static_cast<double>(latticePerSide);

    int particleIndex = 0;
    for (int i = 0; i < latticePerSide && particleIndex < numParticles; i++) {
        for (int j = 0; j < latticePerSide && particleIndex < numParticles; j++) {
            for (int k = 0; k < latticePerSide && particleIndex < numParticles; k++) {
                double x = i * spacing;
                double y = j * spacing;
                double z = k * spacing;
                particles[particleIndex].setPosition(Vector{x, y, z});
                particleIndex++;
            }
        }
    }
    for (int i = 0; i < numBins; i++) {
        radialDistFuncHisto[i] = 0;
    }
    ensembleEnergy = calculateEnergy();
}


// Initialize particles randomly within box of size L
void particleEnsemble::initializeRandom(ntype L) {
    if (L <= 0.0) {
        throw std::invalid_argument("Box length must be positive.");
    }
    boxLength = L;
    for (auto& p : particles) {
        p.random(L);
    }
    ensembleEnergy = calculateEnergy();
}

// Access particle via operator(), with bounds checking
Particle& particleEnsemble::operator()(int i) {
    if (i < 0 || i >= numParticles) {
        throw std::out_of_range("Particle index out of bounds.");
    }
    return particles.at(i);
}

// Initialize a specific particle's position
void particleEnsemble::initParticle(int i, Vector v) {
    if (i < 0 || i >= numParticles) {
        throw std::out_of_range("Particle index out of bounds.");
    }
    particles[i].setPosition(v);
}

// Store current state
void particleEnsemble::store() {
    for (auto& p : particles) {
        p.store();
    }
    pastEnergy = ensembleEnergy;
}

// Restore to stored state
void particleEnsemble::restore() {
    for (auto& p : particles) {
        p.restore();
    }
    ensembleEnergy = pastEnergy;
}

// Get number of particles
int particleEnsemble::getNumParticles() const {
    return numParticles;
}

// Get total energy
ntype particleEnsemble::getEnergy() const {
    return ensembleEnergy;
}

// Set total energy
void particleEnsemble::setEnergy(ntype newEnergy) {
    ensembleEnergy = newEnergy;
}

void particleEnsemble::setOldEnergy(ntype oldEnergy) {
    pastEnergy = oldEnergy;
}

void particleEnsemble::restoreEnergy() {
    ensembleEnergy = pastEnergy;
}   

void particleEnsemble::restoreParticle(int i) {
    particles[i].restore();
}

void particleEnsemble::updateEnsemble(int i, ntype deltaE) {
    ensembleEnergy += deltaE;
    particles[i].store();
    pastEnergy = ensembleEnergy;
}
// Calculate total energy
ntype particleEnsemble::calculateEnergy() {
    ntype energy = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        for (int j = 0; j < i; ++j) { 
            if (i == j) continue;
            Vector r_i = particles[i].getPosition();
            Vector r_j = particles[j].getPosition();
            double distance = potential.computeDistance(r_i, r_j);
            energy += potential.lennardJones(distance);
        }
    }
    ensembleEnergy = energy;
    return energy;
}

ntype particleEnsemble::calculateEnergyDifference(int particleIndex) {
    ntype deltaEnergy = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        if (i == particleIndex) continue;
        Vector r_i = particles[i].getPosition();
        Vector r_jOld = particles[particleIndex].getOldPosition();
        Vector r_j = particles[particleIndex].getPosition();
        double distance = potential.computeDistance(r_i, r_j);
        double oldDistance = potential.computeDistance(r_i, r_jOld);

        deltaEnergy += potential.lennardJones(distance) - potential.lennardJones(oldDistance);
    }
    // std::cout << "Delta diff:" << deltaEnergy << std::endl;
    return deltaEnergy;
}

// Update energy by delta
void particleEnsemble::updateEnergy(ntype deltaEnergy) {
    ensembleEnergy += deltaEnergy; 
}

void particleEnsemble::updateParticle(int i, Vector v) {
    particles[i].setPosition(v);
}

// Get box length
ntype particleEnsemble::getBoxLength() const {
    return boxLength;
}

// Save snapshot of positions to file
void particleEnsemble::saveSnapshotPosition(std::ofstream& file, int step) {
    if (!file.is_open()) {
        std::cerr << "Error: File is not open for writing." << std::endl;
        return;
    }
    file << "Step " << step << std::endl;

    // Write the positions of all particles
    for (auto& p : particles) {
        Vector pos = p.getPosition();
        for (int i = 0; i < 3; ++i) { // Assuming 3 dimensions
            file << pos.get(i) << " ";
        }
        file << std::endl;
    }
}

void particleEnsemble::clearParticles() {
    particles.clear();
    numParticles = 0;
    ensembleEnergy = 0.0;
    pastEnergy = 0.0;
}

void particleEnsemble::addParticle(const Particle& p) {
    particles.push_back(p);
    numParticles++;
}

Vector particleEnsemble::applyMinimumImageConvention(const Vector& position) const {
    Vector newPosition = position;
    for (int i = 0; i < 3; ++i) { 
        if (newPosition.get(i) < 0.0) {
            newPosition.set(i, newPosition.get(i) + boxLength);
        } else if (newPosition.get(i) >= boxLength) {
            newPosition.set(i, newPosition.get(i) - boxLength);
        }
    }
    return newPosition;
}

void particleEnsemble::moveParticle(int i, const Vector& dr) { 
    Vector newPosition = particles[i].getPosition() + dr;
    newPosition = applyMinimumImageConvention(newPosition);
    particles[i].setPosition(newPosition);
}