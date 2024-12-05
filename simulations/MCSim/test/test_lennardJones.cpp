#include <iostream>
#include "particles.hpp"
#include "interactionPotentials.hpp"
#include "vec.hpp"

void testLennardJones() {
    std::cout << "Testing Lennard-Jones potential computation:" << std::endl;

    // Create an ensemble of 100 particles
    int N = 100;
    ntype L = 10.0;
    particleEnsemble ensemble(N, L);

    // Initialize particle positions manually
    for (int i = 0; i < N; ++i) {
        ntype x = static_cast<ntype>(rand()) / RAND_MAX * L;
        ntype y = static_cast<ntype>(rand()) / RAND_MAX * L;
        ntype z = static_cast<ntype>(rand()) / RAND_MAX * L;
        ensemble.initParticle(i, Vector{x, y, z});
    }

    // Manually compute the Lennard-Jones potential
    ntype manualTotal = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            Vector r_i = ensemble(i).getPosition();
            Vector r_j = ensemble(j).getPosition();
            ntype distance = (r_i - r_j).modulus();
            ntype V = 4.0 * (std::pow(1.0 / distance, 12) - std::pow(1.0 / distance, 6));
            manualTotal += V;
        }
    }

    // Compute the Lennard-Jones potential using the function
    interactionPotential potential;
    ntype functionTotal = ensemble.calculateEnergy(); // Assumes calculateEnergy uses lennardJones with sigma and epsilon set to 1

    std::cout << "Manual Lennard-Jones total: " << manualTotal << std::endl;
    std::cout << "Function Lennard-Jones total: " << functionTotal << std::endl;

    if (std::abs(manualTotal - functionTotal) < 1e-5) {
        std::cout << "Lennard-Jones potential test passed." << std::endl;
    } else {
        std::cout << "Lennard-Jones potential test failed." << std::endl;
    }
}

int main() {
    testLennardJones();
    return 0;
}