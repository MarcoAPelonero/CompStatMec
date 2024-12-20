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
    computeVirial();
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
    if (x > boxLength) x -= boxLength;
    if (y < 0) y += boxLength;
    if (y > boxLength) y -= boxLength;
    if (z < 0) z += boxLength;
    if (z > boxLength) z -= boxLength;

    p.setPosition(Vec(x, y, z));
}

Vec ParticleEnsemble::minimumImageConvention(const Particle &p1, const Particle &p2) const {
    Vec r = p1.getPosition() - p2.getPosition();
    double x = r.getX();
    double y = r.getY();
    double z = r.getZ();

    if (x > boxLength / 2) x -= boxLength;
    if (x < -boxLength / 2) x += boxLength;
    if (y > boxLength / 2) y -= boxLength;
    if (y < -boxLength / 2) y += boxLength;
    if (z > boxLength / 2) z -= boxLength;
    if (z < -boxLength / 2) z += boxLength;

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

double ParticleEnsemble::computeTemperature() {
    double temperature = (2.0 / 3.0) * totalKineticEnergy / numParticles;
    return temperature;
}

double ParticleEnsemble::computePressure() {
    double volume = boxLength * boxLength * boxLength;
    double pressure = (2.0 * totalKineticEnergy + totalVirialEnergy) / (3.0 * volume);
    return pressure;
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

void ParticleEnsemble::computeVirial() {
    for (auto& particle : particles) {
        particle.setVirialEnergy(0.0);
    }

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];

            Vec r = minimumImageConvention(p1, p2);
            Vec f = LennardJonesForce(p1, p2);

            double virial = r.dot(f);
            
            p1.setVirialEnergy(p1.getVirialEnergy() + 0.5 * virial);
            p2.setVirialEnergy(p2.getVirialEnergy() + 0.5 * virial);
        }
    }
}

void ParticleEnsemble::computeEnergies() {

    for (Particle& p : particles) {
        p.computeTotalEnergy();
    }
}

void ParticleEnsemble::computeEnsembleEnergies() {
    totalVirialEnergy = 0.0;
    totalKineticEnergy = 0.0;
    totalPotentialEnergy = 0.0;
    totalEnergy = 0.0;

    for (const Particle& p : particles) {
        totalVirialEnergy += p.getVirialEnergy();
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

void ParticleEnsemble::ensembleThermalize(double targetTemperature, double dt, double tauT) {
    double currentTemperature = (2.0 / 3.0) * totalKineticEnergy / numParticles;

    double lambda = std::sqrt(1.0 + (dt / tauT) * (targetTemperature / currentTemperature - 1.0));

    for (Particle& p : particles) {
        Vec velocity = p.getVelocity();
        velocity *= lambda; 
        p.setVelocity(velocity); 
    }
    computeKinetics();
    computeEnergies();
    computeEnsembleEnergies();
}

void ParticleEnsemble::ensemblePressurize(double targetPressure, double dt, double tauP) {

    double volume = std::pow(boxLength, 3);
    double currentPressure = (2.0 * totalKineticEnergy + totalVirialEnergy) / (3.0 * volume);

    double lambda = 1.0 - (dt / tauP) * ((targetPressure - currentPressure));

    const double lambda_min = 0.9;
    const double lambda_max = 1.1;
    if (lambda < lambda_min) lambda = lambda_min;
    if (lambda > lambda_max) lambda = lambda_max;

    double boxScale = std::pow(lambda, 1.0 / 3.0);
    boxLength *= boxScale;

     for (auto &particle : particles) {
         Vec position = particle.getPosition();
        position *= boxScale;
        particle.setPosition(position);
        periodicBoundaryConditions(particle);   
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
    computeVirial();
    computeEnergies();
    computeEnsembleEnergies();
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
    computeVirial();
    computeEnergies();
    computeEnsembleEnergies();
}

void ParticleEnsemble::ensembleStepSpeedVerlet(double dt) {
    
    computeForces();

    std::vector<Vec> oldAccelerations;
    oldAccelerations.clear();
    
    for (Particle& p : particles) {
        Vec force = p.getForce();
        Vec velocity = p.getVelocity();
        Vec position = p.getPosition();

        Vec acceleration = force / p.getMass();
        oldAccelerations.push_back(acceleration);

        position += velocity * dt + 0.5 * acceleration * dt * dt;
        p.setPosition(position);

        periodicBoundaryConditions(p);
    }

    computeForces();
    int accelerationCounter = 0;

    for (Particle& p : particles) {
        Vec force = p.getForce();
        Vec acceleration = force / p.getMass();
        Vec velocity = p.getVelocity();
        velocity += (acceleration + oldAccelerations[accelerationCounter]) * 0.5 * dt;
        p.setVelocity(velocity);
        accelerationCounter++;
    }

    computeVirial();
    computePotentials();
    computeKinetics();
    computeEnergies();
}

void ParticleEnsemble::ensembleSnapshot(std::ofstream &file) {
    file << boxLength << "\n";
    for (const auto &particle : particles) {
        file << particle.getPosition() << " "
             << particle.getVelocity() << " "
             << particle.getPotentialEnergy() << " "
             << particle.getKineticEnergy() << " "
             << particle.getTotalEnergy() << "\n";
    }
    file << "\n"; 
}


void ParticleEnsemble::resetRDFHistogram() {
    std::fill(rdfHistogram.begin(), rdfHistogram.end(), 0.0);
}

void ParticleEnsemble::computeRadialDistributionFunctionDirect() {
    rdfCutoff = boxLength / 2.0;
    rdfBinWidth = rdfCutoff / rdfBins;

    resetRDFHistogram();

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Vec dr = minimumImageConvention(particles[i], particles[j]);
            double r = dr.length();
            if (r < rdfCutoff) {
                int bin = (int)std::floor(r / rdfBinWidth);
                if (bin >= 0 && bin < rdfBins) {
                    rdfHistogram[bin] += 2.0; 
                }
            }
        }
    }

    double volume = boxLength * boxLength * boxLength;
    double density = (double)particles.size() / volume;

    for (int i = 0; i < rdfBins; ++i) {
        double rLower = i * rdfBinWidth;
        double rUpper = rLower + rdfBinWidth;
        double shellVolume = (4.0 / 3.0) * M_PI * (rUpper*rUpper*rUpper - rLower*rLower*rLower);
        double idealCount = density * shellVolume * (double)particles.size();
        rdfHistogram[i] /= idealCount;
    }
}

void ParticleEnsemble::buildCells(std::vector<Cell> &cells, int &numCellsPerDim, double &cellSize) const {
    cellSize = rdfCutoff; 
    numCellsPerDim = (int)std::floor(boxLength / cellSize);
    if (numCellsPerDim < 1) numCellsPerDim = 1; 

    cells.clear();
    cells.resize(numCellsPerDim * numCellsPerDim * numCellsPerDim);

    for (size_t i = 0; i < particles.size(); ++i) {
        Vec pos = particles[i].getPosition();

        int cx = (int)std::floor(pos.getX() / cellSize);
        int cy = (int)std::floor(pos.getY() / cellSize);
        int cz = (int)std::floor(pos.getZ() / cellSize);

        if (cx == numCellsPerDim) cx = numCellsPerDim - 1;
        if (cy == numCellsPerDim) cy = numCellsPerDim - 1;
        if (cz == numCellsPerDim) cz = numCellsPerDim - 1;

        int cellIndex = cx + cy * numCellsPerDim + cz * numCellsPerDim * numCellsPerDim;
        cells[cellIndex].particleIndices.push_back((int)i);
    }
}

void ParticleEnsemble::computeRadialDistributionFunctionCellMethod() {
    rdfCutoff = boxLength / 2.0;
    rdfBinWidth = rdfCutoff / rdfBins;

    resetRDFHistogram();

    std::vector<Cell> cells;
    int numCellsPerDim;
    double cellSize;
    buildCells(cells, numCellsPerDim, cellSize);

    int neighborOffsets[3] = {-1, 0, 1};

    for (int cx = 0; cx < numCellsPerDim; ++cx) {
        for (int cy = 0; cy < numCellsPerDim; ++cy) {
            for (int cz = 0; cz < numCellsPerDim; ++cz) {
                int cellIndex = cx + cy * numCellsPerDim + cz * numCellsPerDim * numCellsPerDim;
                const std::vector<int> &cParticles = cells[cellIndex].particleIndices;

                for (int nx : neighborOffsets) {
                    for (int ny : neighborOffsets) {
                        for (int nz : neighborOffsets) {
                            int ncx = (cx + nx + numCellsPerDim) % numCellsPerDim;
                            int ncy = (cy + ny + numCellsPerDim) % numCellsPerDim;
                            int ncz = (cz + nz + numCellsPerDim) % numCellsPerDim;

                            int neighborCellIndex = ncx + ncy * numCellsPerDim + ncz * numCellsPerDim * numCellsPerDim;
                            const std::vector<int> &nParticles = cells[neighborCellIndex].particleIndices;

                            if (neighborCellIndex == cellIndex) {
                                for (size_t i = 0; i < cParticles.size(); ++i) {
                                    for (size_t j = i + 1; j < cParticles.size(); ++j) {
                                        Vec dr = minimumImageConvention(particles[cParticles[i]], particles[cParticles[j]]);
                                        double r = dr.length();
                                        if (r < rdfCutoff) {
                                            int bin = (int)std::floor(r / rdfBinWidth);
                                            if (bin >= 0 && bin < rdfBins) {
                                                rdfHistogram[bin] += 2.0; 
                                            }
                                        }
                                    }
                                }
                            } else {
                                // Different cell
                                for (auto pi : cParticles) {
                                    for (auto pj : nParticles) {
                                        if (pi == pj) continue; 
                                        Vec dr = minimumImageConvention(particles[pi], particles[pj]);
                                        double r = dr.length();
                                        if (r < rdfCutoff) {
                                            int bin = (int)std::floor(r / rdfBinWidth);
                                            if (bin >= 0 && bin < rdfBins) {
                                                rdfHistogram[bin] += 1.0;
                                            }
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    // Normalize
    double volume = boxLength * boxLength * boxLength;
    double density = (double)particles.size() / volume;

    for (int i = 0; i < rdfBins; ++i) {
        double rLower = i * rdfBinWidth;
        double rUpper = rLower + rdfBinWidth;
        double shellVolume = (4.0 / 3.0) * M_PI * (rUpper*rUpper*rUpper - rLower*rLower*rLower);
        double idealCount = density * shellVolume * (double)particles.size();
        rdfHistogram[i] /= idealCount;
    }
}

void ParticleEnsemble::computeThermodynamicsFromRDF(std::ofstream &file) {
    double temperature = computeTemperature();
    double volume = boxLength * boxLength * boxLength;
    double density = static_cast<double>(numParticles) / volume;

    double dr = rdfBinWidth;
    double energy_per_particle = 0.0;
    double pressure_term = 0.0;
    double k_B = 1.0; 

    for (int i = 0; i < rdfBins; ++i) {
        double r = (i + 0.5) * dr;   
        double g_r = rdfHistogram[i];

        double r2 = r*r;
        double inv_r2 = 1.0 / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;   
        double inv_r12 = inv_r6 * inv_r6;           

        double u_r = 4.0 * (inv_r12 - inv_r6);     

        energy_per_particle += g_r * u_r * (r2) * dr;

        double inv_r7 = inv_r6 / (r*r);  
        double inv_r13 = inv_r12 / (r*r); 
        double du_dr = 24.0 * (2.0 * inv_r13 - inv_r7);

        pressure_term += g_r * (4.0 * M_PI * (r*r*r) * du_dr) * dr;
    }

    energy_per_particle *= 2.0 * M_PI * density;
    double total_potential_energy = energy_per_particle * numParticles;

    double p = density * k_B * temperature - (2.0 / 3.0) * M_PI * density * density * pressure_term;

    file << temperature << " " << energy_per_particle << " " << total_potential_energy << " " << p << "\n";
}


void ParticleEnsemble::printRadialDistributionFunction(std::ofstream &file) {
    for (int i = 0; i < rdfBins; ++i) {
        double r = (i + 0.5) * rdfBinWidth;
        file << r << " " << rdfHistogram[i] << "\n";
    }
    file << "\n";
}