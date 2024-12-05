#include <iostream>
#include "particles.hpp"
#include "configuration.hpp"
#include "interactionPotentials.hpp"
#include "vec.hpp"

void testTwoIndividualParticles() {
    std::cout << "Testing two individual particles:" << std::endl;

    // Create two particles with specific positions
    Vector pos1{0.0, 0.0, 0.0};
    Vector pos2{1.0, 0.0, 0.0}; // 1 unit apart

    Particle p1(pos1);
    Particle p2(pos2);

    // Calculate the distance between the particles
    ntype distance = (p1.getPosition() - (p2.getPosition())).modulus();
    std::cout << "Distance between particles: " << distance << std::endl;

    // Calculate the Lennard-Jones potential
    interactionPotential potential;
    ntype ljPotential = potential.lennardJones(distance);
    std::cout << "Lennard-Jones potential: " << ljPotential << std::endl;

    // Verify the potential is as expected
    ntype expectedPotential = 4.0 * 1.0 * (std::pow(0.9 / distance, 12) - std::pow(0.9 / distance, 6));
    if (std::abs(ljPotential - expectedPotential) < 1e-5) {
        std::cout << "Lennard-Jones potential test passed." << std::endl;
    } else {
        std::cout << "Lennard-Jones potential test failed." << std::endl;
    }
}

void testTwoParticlesInEnsemble() {
    std::cout << "Testing two particles in an ensemble:" << std::endl;

    // Create two particles with specific positions
    Vector pos1{0.0, 0.0, 0.0};
    Vector pos2{1.57345, 0.0, 0.0}; // 1 unit apart

    // Create an ensemble and add the particles
    particleEnsemble ensemble(2, 10.0);
    ensemble.initParticle(0, pos1);
    ensemble.initParticle(1, pos2);

    // Calculate the distance between the particles
    ntype distance = ensemble.computeDistance(ensemble(0).getPosition(), ensemble(1).getPosition());
    std::cout << "Distance between particles: " << distance << std::endl;

    // Calculate the Lennard-Jones potential
    interactionPotential potential;
    ntype ljPotential = potential.lennardJones(distance);
    std::cout << "Lennard-Jones potential: " << ljPotential << std::endl;

    // Verify the potential is as expected
    ntype expectedPotential = 4.0 * 1.0 * (std::pow(0.9 / distance, 12) - std::pow(0.9 / distance, 6));
    std::cout << "Distance between particles: " << distance << std::endl;
    std::cout << "Lennard-Jones potential: " << ljPotential << std::endl;
    std::cout << "Expected potential: " << expectedPotential << std::endl;
    if (std::abs(ljPotential - expectedPotential) < 1e-5) {
        std::cout << "Lennard-Jones potential test passed." << std::endl;
    } else {
        std::cout << "Lennard-Jones potential test failed." << std::endl;
    }
}

int main() {
    testTwoIndividualParticles();
    testTwoParticlesInEnsemble();
    return 0;
}