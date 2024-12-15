#include "particleEnsemble.hpp" // Include the header file

void ParticleEnsemble::initializeParticles() {
    int numPerSide = std::ceil(std::cbrt(numParticles)); 
    double spacing = boxLength / numPerSide; 

    int particleCount = 0;
    for (int x = 0; x < numPerSide && particleCount < numParticles; ++x) {
        for (int y = 0; y < numPerSide && particleCount < numParticles; ++y) {
            for (int z = 0; z < numPerSide && particleCount < numParticles; ++z) {
                Vec position(x * spacing, y * spacing, z * spacing);
                Vec velocity(gaussian(rng), gaussian(rng), gaussian(rng));
                Particle particle(1.0, position, velocity); 
                particles.push_back(particle);
                ++particleCount;
            }
        }
    }

    computePotentials();
    computeKinetics();
    computeEnergies();
    computeEnsembleEnergies();
}

std::vector<Particle>& ParticleEnsemble::getParticles() {
    return particles;
}

const std::vector<Particle>& ParticleEnsemble::getParticles() const {
    return particles;
}

void ParticleEnsemble::showParticles() const {
    for (size_t i = 0; i < particles.size(); ++i) {
        std::cout << "Particle " << i + 1 << ":\n";
        particles[i].show("\t");
    }
}

void ParticleEnsemble::periodicBoundaryConditions(Particle &p) const {
    Vec pos = p.getPosition();
    double x = pos.getX();
    double y = pos.getY();
    double z = pos.getZ();

    if (x < 0) x += boxLength;
    if (x >= boxLength) x -= boxLength;
    if (y < 0) y += boxLength;
    if (y >= boxLength) y -= boxLength;
    if (z < 0) z += boxLength;
    if (z >= boxLength) z -= boxLength;

    p.setPosition(Vec(x, y, z));
}

Vec ParticleEnsemble::minimumImageConvention(const Particle &p1, const Particle &p2) const {
    Vec r = p1.getPosition() - p2.getPosition();
    double x = r.getX();
    double y = r.getY();
    double z = r.getZ();

    if (x > boxLength / 2) x -= boxLength;
    if (x <= -boxLength / 2) x += boxLength;
    if (y > boxLength / 2) y -= boxLength;
    if (y <= -boxLength / 2) y += boxLength;
    if (z > boxLength / 2) z -= boxLength;
    if (z <= -boxLength / 2) z += boxLength;

    return Vec(x, y, z);
}

double ParticleEnsemble::LennardJonesPotential(const Particle &p1, const Particle &p2) const {
    Vec r = minimumImageConvention(p1, p2); 
    double r_len = r.length();
    
    if (r_len == 0.0 || r_len > 4.0) return 0.0; 
    double r2 = r_len * r_len;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    return 4 * (1 / r12 - 1 / r6);
}

void ParticleEnsemble::computePotentials() {
    for (auto& particle : particles) {
        particle.setPotentialEnergy(0.0);
    }

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];
            double potential = LennardJonesPotential(p1, p2);
            
            p1.setPotentialEnergy(p1.getPotentialEnergy() + 0.5 * potential);
            p2.setPotentialEnergy(p2.getPotentialEnergy() + 0.5 * potential);
        }
    }
}

void ParticleEnsemble::computeKinetics() {
    for (Particle& p : particles) {
        p.computeKineticEnergy();
    }
}

void ParticleEnsemble::computeEnergies() {

    for (Particle& p : particles) {
        p.computeTotalEnergy();
    }
}

void ParticleEnsemble::computeEnsembleEnergies() {
    totalKineticEnergy = 0.0;
    totalPotentialEnergy = 0.0;
    totalEnergy = 0.0;

    for (const Particle& p : particles) {
        totalKineticEnergy += p.getKineticEnergy();
        totalPotentialEnergy += p.getPotentialEnergy();
        totalEnergy += p.getTotalEnergy();
    }
}

Vec ParticleEnsemble::LennardJonesForce(const Particle &p1, const Particle &p2) const {
    Vec r = minimumImageConvention(p1, p2); 
    double r_len = r.length();
    
    if (r_len == 0.0) return Vec(0.0, 0.0, 0.0);
    if (r_len > 4.0) return Vec(0.0, 0.0, 0.0); 

    double r2 = r_len * r_len;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    double force_magnitude = 24 * (2 / r12 - 1 / r6) / r2;

    return force_magnitude * r; 
}


void ParticleEnsemble::computeForces() {
    for (Particle& p : particles) {
        p.setForce(Vec(0, 0, 0));
    }

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];
            Vec force = LennardJonesForce(p1, p2);
            p1.addForce(force);
            p2.addForce(-force);
        }
    }
}

void ParticleEnsemble::ensembleThermalize(double targetTemperature) {
    computeKinetics();

    double currentTemperature = (2.0 / 3.0) * totalKineticEnergy / numParticles;

    double scale = std::sqrt(targetTemperature / currentTemperature);

    for (Particle& p : particles) {
        Vec velocity = p.getVelocity();
        velocity *= scale; 
        p.setVelocity(velocity); 
    }

    computeKinetics();
    computeEnergies();
    computeEnsembleEnergies();
}

void ParticleEnsemble::ensembleStepEuler(double dt) {
    
    computeForces();

    for (Particle& p : particles) {
        Vec force = p.getForce();
        Vec acceleration = force / p.getMass();
        Vec velocity = p.getVelocity();
        Vec position = p.getPosition();
        position += velocity * dt;
        velocity += acceleration * dt;
        p.setVelocity(velocity);
        p.setPosition(position);

        periodicBoundaryConditions(p);
    }

    computePotentials();
    computeKinetics();
    computeEnergies();
}

void ParticleEnsemble::ensembleStepEulerCromer(double dt) {
    
    computeForces();
    
    for (Particle& p : particles) {
        Vec force = p.getForce();
        Vec acceleration = force / p.getMass();
        Vec velocity = p.getVelocity();
        Vec position = p.getPosition();
        velocity += acceleration * dt;
        position += velocity * dt;
        p.setVelocity(velocity);
        p.setPosition(position);

        periodicBoundaryConditions(p);
    }

    computePotentials();
    computeKinetics();
    computeEnergies();
}

void ParticleEnsemble::ensembleSnapshot(std::ofstream &file) {
    for (const auto &particle : particles) {
        file << particle.getPosition() << " "
             << particle.getVelocity() << " "
             << particle.getPotentialEnergy() << " "
             << particle.getKineticEnergy() << " "
             << particle.getTotalEnergy() << "\n";
    }
    file << "\n"; // Skip a line after each snapshot
}