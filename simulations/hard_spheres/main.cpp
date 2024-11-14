// main.cpp

#include <iostream>
#include <cmath>
#include <random>
#include <initializer_list>
#include <string>
#include <fstream> // Include for file operations
#include "vector.hpp"

constexpr int N_PARTICLES = 150;        // Number of particles
constexpr double BOX_SIZE = 10.0;       // Size of the simulation box
constexpr double DIAMETER = 2.0;        // Diameter of hard spheres
constexpr double INFINITE_ENERGY = 1e10; // Representation of infinite energy
constexpr int N_STEPS = 100000;         // Number of simulation steps
constexpr double MAX_DISPLACEMENT = 0.1; // Maximum displacement per step
constexpr double BOLTZMANN_CONSTANT = 1.0; // Boltzmann constant
constexpr double TEMPERATURE = 1.0;     // Temperature of the system
constexpr int OUTPUT_INTERVAL = 1000;   // Interval for outputting positions
constexpr int DIM = 3;                  // Dimension of the system


// Potential base class
class Potential {
public:
    virtual double compute(Vector& r1, Vector& r2) const = 0;
    virtual ~Potential() {}
};

// Hard Sphere Potential class
class HardSpherePotential : public Potential {
public:
    double compute(Vector& r1, Vector& r2) const override {
        double distance = (r1 - r2).modulus();
        return (distance < DIAMETER) ? INFINITE_ENERGY : 0.0;
    }
};

// Function to apply periodic boundary conditions
void applyPeriodicBoundary(Vector& position, double boxSize) {
    for (int i = 0; i < DIM; ++i) {
        if (position(i) >= boxSize) {
            position(i) -= boxSize;
        } else if (position(i) < 0.0) {
            position(i) += boxSize;
        }
    }
}

// Simulation function
void runSimulation(Potential& potential, Vector positions[N_PARTICLES], std::ofstream& outputFile) {
    int accepted = 0;
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::uniform_int_distribution<int> particleDist(0, N_PARTICLES - 1);

    // Record initial positions
    for (int i = 0; i < N_PARTICLES; ++i) {
        for (int d = 0; d < DIM; ++d) {
            outputFile << positions[i](d) << " ";
        }
        outputFile << "\n";
    }
    outputFile << "\n"; // Separate time steps

    for (int step = 0; step < N_STEPS; ++step) {
        // Select a particle at random
        int i = particleDist(rng);

        // Propose a new position
        Vector displacement;
        for (int d = 0; d < DIM; ++d) {
            std::uniform_real_distribution<double> dispDist(-MAX_DISPLACEMENT, MAX_DISPLACEMENT);
            displacement(d) = dispDist(rng);
        }
        Vector newPosition = positions[i] + displacement;

        // Apply periodic boundary conditions
        applyPeriodicBoundary(newPosition, BOX_SIZE);

        // Compute energy difference
        double deltaE = 0.0;
        for (int j = 0; j < N_PARTICLES; ++j) {
            if (j != i) {
                double oldEnergy = potential.compute(positions[i], positions[j]);
                double newEnergy = potential.compute(newPosition, positions[j]);
                deltaE += newEnergy - oldEnergy;
            }
        }

        // Metropolis criterion
        if (deltaE <= 0.0 || exp(-deltaE / (BOLTZMANN_CONSTANT * TEMPERATURE)) > dist(rng)) {
            positions[i] = newPosition;
            accepted++;
        }

        // Output positions at intervals
        if ((step + 1) % OUTPUT_INTERVAL == 0) {
            for (int i = 0; i < N_PARTICLES; ++i) {
                for (int d = 0; d < DIM; ++d) {
                    outputFile << positions[i](d) << " ";
                }
                outputFile << "\n";
            }
            outputFile << "\n"; // Separate time steps
        }
    }
}

// Main function
int main() {
    // Initialize particle positions
    Vector positions[N_PARTICLES];

    // Random number generator for initial positions
    std::mt19937 rng;
    std::uniform_real_distribution<double> positionDist(0.0, BOX_SIZE);

    for (int i = 0; i < N_PARTICLES; ++i) {
        for (int d = 0; d < DIM; ++d) {
            positions[i](d) = positionDist(rng);
        }
    }

    // Instantiate the hard sphere potential
    HardSpherePotential potential;

    // Open output file
    std::ofstream outputFile("positions.dat");

    // Run the simulation
    runSimulation(potential, positions, outputFile);

    // Close the output file
    outputFile.close();

    return 0;
}
