// lennard_jones_sim.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <cstdlib>

// --- Basic 3D vector class ---
class Vector {
public:
    double x, y, z;
    Vector(double x_=0.0, double y_=0.0, double z_=0.0) : x(x_), y(y_), z(z_) {}

    Vector operator+(const Vector &other) const {
        return Vector(x + other.x, y + other.y, z + other.z);
    }
    Vector operator-(const Vector &other) const {
        return Vector(x - other.x, y - other.y, z - other.z);
    }
    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }
};

// --- Particle class ---
class Particle {
public:
    Vector pos;
    Particle(const Vector &p) : pos(p) {}
};

// --- Periodic Boundary Conditions with Minimum Image Convention ---
double minImage(double d, double L) {
    return d - L * std::round(d / L);
}

double distance(const Vector &a, const Vector &b, double L) {
    double dx = minImage(a.x - b.x, L);
    double dy = minImage(a.y - b.y, L);
    double dz = minImage(a.z - b.z, L);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// --- Lennard-Jones potential (epsilon = 1, sigma = 1) ---
double ljPotential(double r) {
    if(r == 0) return 1e12; // avoid singularity
    double r2 = r * r;
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;
    return 4.0 * (inv_r12 - inv_r6);
}

// --- Full potential energy difference ---
double computeEnergyDiffFull(const std::vector<Particle> &particles, int i, const Vector &newPos, double L) {
    double dE = 0.0;
    for (int j = 0; j < particles.size(); j++) {
        if (j == i) continue;
        double r_old = distance(particles[i].pos, particles[j].pos, L);
        double r_new = distance(newPos, particles[j].pos, L);
        dE += ljPotential(r_new) - ljPotential(r_old);
    }
    return dE;
}

// --- Truncated potential energy difference ---
double computeEnergyDiffTruncated(const std::vector<Particle> &particles, int i, const Vector &newPos, double L, double cutoff) {
    double dE = 0.0;
    for (int j = 0; j < particles.size(); j++) {
        if (j == i) continue;
        double r_old = distance(particles[i].pos, particles[j].pos, L);
        double r_new = distance(newPos, particles[j].pos, L);
        double pot_old = (r_old < cutoff) ? ljPotential(r_old) : 0.0;
        double pot_new = (r_new < cutoff) ? ljPotential(r_new) : 0.0;
        dE += pot_new - pot_old;
    }
    return dE;
}

// --- Helper for linked cells: build cell list ---
void buildCellList(const std::vector<Particle> &particles, double L, double cutoff,
                   std::vector< std::vector<int> > &cellList, int &nCells) {
    nCells = static_cast<int>(L / cutoff);
    if (nCells < 1) nCells = 1;
    int totalCells = nCells * nCells * nCells;
    cellList.clear();
    cellList.resize(totalCells);
    for (int i = 0; i < particles.size(); i++) {
        int ix = static_cast<int>(particles[i].pos.x / L * nCells) % nCells;
        int iy = static_cast<int>(particles[i].pos.y / L * nCells) % nCells;
        int iz = static_cast<int>(particles[i].pos.z / L * nCells) % nCells;
        int cellIndex = ix + iy * nCells + iz * nCells * nCells;
        cellList[cellIndex].push_back(i);
    }
}

// --- Linked cells energy difference ---
double computeEnergyDiffLinked(const std::vector<Particle> &particles, int i, const Vector &newPos, double L, double cutoff) {
    int nCells;
    std::vector< std::vector<int> > cellList;
    buildCellList(particles, L, cutoff, cellList, nCells);

    // Old energy using neighboring cells of particle i's current cell
    int ix_old = static_cast<int>(particles[i].pos.x / L * nCells) % nCells;
    int iy_old = static_cast<int>(particles[i].pos.y / L * nCells) % nCells;
    int iz_old = static_cast<int>(particles[i].pos.z / L * nCells) % nCells;
    double energy_old = 0.0;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                int ix = (ix_old + dx + nCells) % nCells;
                int iy = (iy_old + dy + nCells) % nCells;
                int iz = (iz_old + dz + nCells) % nCells;
                int cellIndex = ix + iy * nCells + iz * nCells * nCells;
                for (int j : cellList[cellIndex]) {
                    if (j == i) continue;
                    double r = distance(particles[i].pos, particles[j].pos, L);
                    if (r < cutoff)
                        energy_old += ljPotential(r);
                }
            }
        }
    }

    // New energy using the candidate new position
    int ix_new = static_cast<int>(newPos.x / L * nCells) % nCells;
    int iy_new = static_cast<int>(newPos.y / L * nCells) % nCells;
    int iz_new = static_cast<int>(newPos.z / L * nCells) % nCells;
    double energy_new = 0.0;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                int ix = (ix_new + dx + nCells) % nCells;
                int iy = (iy_new + dy + nCells) % nCells;
                int iz = (iz_new + dz + nCells) % nCells;
                int cellIndex = ix + iy * nCells + iz * nCells * nCells;
                for (int j : cellList[cellIndex]) {
                    if (j == i) continue;
                    double r = distance(newPos, particles[j].pos, L);
                    if (r < cutoff)
                        energy_new += ljPotential(r);
                }
            }
        }
    }
    return energy_new - energy_old;
}

// --- Neighbor list energy difference ---
double computeEnergyDiffNeighbor(const std::vector<Particle> &particles, int i, const Vector &newPos, double L, double cutoff) {
    std::vector<int> neighbors;
    for (int j = 0; j < particles.size(); j++) {
        if (j == i) continue;
        double r = distance(particles[i].pos, particles[j].pos, L);
        if (r < cutoff) {
            neighbors.push_back(j);
        }
    }
    double energy_old = 0.0, energy_new = 0.0;
    for (int j : neighbors) {
        double r_old = distance(particles[i].pos, particles[j].pos, L);
        double r_new = distance(newPos, particles[j].pos, L);
        energy_old += ljPotential(r_old);
        energy_new += ljPotential(r_new);
    }
    return energy_new - energy_old;
}

// --- Apply periodic boundary conditions (using minimum image convention) ---
Vector applyPBC(const Vector &v, double L) {
    double x = v.x - L * std::floor(v.x / L);
    double y = v.y - L * std::floor(v.y / L);
    double z = v.z - L * std::floor(v.z / L);
    return Vector(x, y, z);
}

// --- Main simulation ---
int main(int argc, char* argv[]) {
    // Usage: ./lennard_jones_sim mode N steps [seed]
    if(argc < 4) {
        std::cerr << "Usage: " << argv[0] << " mode N steps [seed]\n";
        std::cerr << " mode: 0=full, 1=truncated, 2=linked cells, 3=neighbor list\n";
        return 1;
    }
    int mode = std::atoi(argv[1]);
    int N = std::atoi(argv[2]);
    int steps = std::atoi(argv[3]);
    unsigned int seed = (argc >= 5) ? std::atoi(argv[4]) : 12345;
    
    // Simulation parameters
    double T = 1.0;             // Temperature (reduced units: kT = 1)
    double delta = 0.1;         // Maximum displacement
    double cutoff = 2.5;        // Cutoff distance for truncated, linked, and neighbor list modes
    double density = 0.7;       // Fixed density as requested
    double L = std::cbrt(N / density);  // Box length computed from N and density

    // Setup random number generators
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    std::uniform_real_distribution<double> disp(-delta, delta);

    // Initialize particles with random positions inside the box
    std::vector<Particle> particles;
    particles.reserve(N);
    for (int i = 0; i < N; i++) {
        double x = uniform01(rng) * L;
        double y = uniform01(rng) * L;
        double z = uniform01(rng) * L;
        particles.push_back(Particle(Vector(x, y, z)));
    }

    // Start simulation timer
    auto start_time = std::chrono::high_resolution_clock::now();

    int accepted = 0;
    // Monte Carlo Metropolis loop
    for (int step = 0; step < steps; step++) {
        int i = rng() % N; // select random particle
        Vector oldPos = particles[i].pos;
        Vector move(disp(rng), disp(rng), disp(rng));
        Vector newPos = applyPBC(oldPos + move, L);

        // Compute energy difference based on simulation mode
        double dE = 0.0;
        switch(mode) {
            case 0:
                dE = computeEnergyDiffFull(particles, i, newPos, L);
                break;
            case 1:
                dE = computeEnergyDiffTruncated(particles, i, newPos, L, cutoff);
                break;
            case 2:
                dE = computeEnergyDiffLinked(particles, i, newPos, L, cutoff);
                break;
            case 3:
                dE = computeEnergyDiffNeighbor(particles, i, newPos, L, cutoff);
                break;
            default:
                std::cerr << "Invalid mode.\n";
                return 1;
        }
        // Metropolis acceptance criterion
        if(dE <= 0 || std::exp(-dE / T) > uniform01(rng)) {
            particles[i].pos = newPos;
            accepted++;
        }
    }

    // End simulation timer
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    double avgTime = elapsed.count() / steps;

    // Output: mode,N,steps,accepted,avg_time_per_step
    std::cout << mode << "," << N << "," << steps << "," << accepted << "," << avgTime << "\n";
    return 0;
}
