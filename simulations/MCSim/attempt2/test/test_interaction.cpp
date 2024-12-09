#include <iostream>
#include <cmath>
#include <vector>
#include "interactionPotential.hpp"
#include "particle.hpp"

void testDefaultConstructor() {
    interactionPotential ip;
    if (ip.getSigma() == 1.0 && ip.getEpsilon() == 1.0 && ip.getL() == 10.0 && ip.getRcut() == 2.5) {
        std::cout << "Default constructor test passed." << std::endl;
    } else {
        std::cerr << "Default constructor test failed." << std::endl;
    }
}

void testParameterizedConstructor() {
    ntype boxLength = 10.0;
    interactionPotential ip(boxLength);
    if (ip.getL() == boxLength && ip.getSigma() == 1.0 && ip.getEpsilon() == 1.0 && ip.getRcut() == 2.5) {
        std::cout << "Parameterized constructor test passed." << std::endl;
    } else {
        std::cerr << "Parameterized constructor test failed." << std::endl;
    }
}

void testFullConstructor() {
    ntype sigma = 1.0;
    ntype epsilon = 1.0;
    ntype boxLength = 10.0;
    ntype rcut = 2.5;
    interactionPotential ip(sigma, epsilon, boxLength, rcut);
    if (ip.getSigma() == sigma && ip.getEpsilon() == epsilon && ip.getL() == boxLength && ip.getRcut() == rcut) {
        std::cout << "Full constructor test passed." << std::endl;
    } else {
        std::cerr << "Full constructor test failed." << std::endl;
    }
}

void testComputeDistance() {
    interactionPotential ip;
    Vector r1{1.0, 2.0, 3.0};
    Vector r2{4.0, 5.0, 6.0};
    ntype distance = ip.computeDistance(r1, r2);
    if (std::abs(distance - std::sqrt(27.0)) < 1e-6) {
        std::cout << "Compute distance test passed." << std::endl;
    } else {
        std::cerr << "Compute distance test failed." << std::endl;
    }
}

void testLennardJones() {
    interactionPotential ip;
    ntype r = 1.0;
    ntype lj = ip.lennardJones(r);
    ntype expectedLJ = 4 * (std::pow(1.0 / r, 12) - std::pow(1.0 / r, 6));
    if (std::abs(lj - expectedLJ) < 1e-6) {
        std::cout << "Lennard-Jones potential test passed." << std::endl;
    } else {
        std::cerr << "Lennard-Jones potential test failed." << std::endl;
    }
}

void testCutLennardJones() {
    interactionPotential ip;
    ntype r = 1.0;
    ntype clj = ip.cutLennardJones(r);
    ntype expectedCLJ = ip.lennardJones(r) - ip.lennardJones(ip.getRcut());
    if (std::abs(clj - expectedCLJ) < 1e-6) {
        std::cout << "Cut Lennard-Jones potential test passed." << std::endl;
    } else {
        std::cerr << "Cut Lennard-Jones potential test failed." << std::endl;
    }
}

void testInteractionWithParticles() {
    ntype sigma = 1.0;
    ntype epsilon = 1.0;
    ntype boxLength = 10.0;
    ntype rcut = 2.5;
    interactionPotential ip(sigma, epsilon, boxLength, rcut);

    Particle p1(Vector{1.0, 2.0, 3.0});
    Particle p2(Vector{4.0, 5.0, 6.0});
    Particle p3(Vector{7.0, 8.0, 9.0});
    Particle p4(Vector{10.0, 11.0, 12.0});

    std::vector<Particle> particles = {p1, p2, p3, p4};
    
    bool allTestsPassed = true;

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            ntype distance = ip.computeDistance(particles[i].getPosition(), particles[j].getPosition());
            
            // Test for regular Lennard-Jones potential
            ntype lj = ip.lennardJones(distance);
            ntype expectedLJ = 4 * (std::pow(sigma / distance, 12) - std::pow(sigma / distance, 6));

            if (std::abs(lj - expectedLJ) >= 1e-6) {
                allTestsPassed = false;
                std::cerr << "Lennard-Jones potential test failed between particle " 
                          << i << " and particle " << j << std::endl;
            }

            ntype clj = ip.cutLennardJones(distance);

            if (distance >= rcut) {
                if (std::abs(clj) > 1e-6) {
                    allTestsPassed = false;
                    std::cerr << "Cut Lennard-Jones potential test failed between particle " 
                              << i << " and particle " << j << " (distance: " << distance << ")" << std::endl;
                }
            } else {
                ntype expectedCLJ = 4 * (std::pow(1.0 / distance, 12) - std::pow(1.0 / distance, 6));
                if (std::abs(clj - expectedCLJ) >= 1e-6) {
                    allTestsPassed = false;
                    std::cerr << "Cut Lennard-Jones potential test failed between particle " 
                              << i << " and particle " << j << std::endl;
                }
            }
        }
    }

    if (allTestsPassed) {
        std::cout << "Multiple particles interaction test passed." << std::endl;
    } else {
        std::cerr << "Multiple particles interaction test failed." << std::endl;
    }
}


void testMultipleParticles() {
    ntype sigma = 1.0;
    ntype epsilon = 1.0;
    ntype boxLength = 10.0;
    ntype rcut = 2.5;
    interactionPotential ip(sigma, epsilon, boxLength, rcut);

    std::vector<Particle> particles;
    particles.emplace_back(Vector{1.0, 2.0, 3.0});
    particles.emplace_back(Vector{4.0, 5.0, 6.0});
    particles.emplace_back(Vector{7.0, 8.0, 9.0});
    particles.emplace_back(Vector{10.0, 11.0, 12.0});

    bool allTestsPassed = true;
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            ntype distance = ip.computeDistance(particles[i].getPosition(), particles[j].getPosition());
            ntype lj = ip.lennardJones(distance);
            ntype clj = ip.cutLennardJones(distance);

            ntype expectedLJ = 4 * (std::pow(1.0 / distance, 12) - std::pow(1.0 / distance, 6));
            ntype expectedCLJ = ip.lennardJones(distance) - ip.lennardJones(ip.getRcut());

            if (std::abs(lj - expectedLJ) >= 1e-6 || std::abs(clj - expectedCLJ) >= 1e-6) {
                allTestsPassed = false;
                std::cerr << "Test failed between particle " << i << " and particle " << j << std::endl;
            }
        }
    }

    if (allTestsPassed) {
        std::cout << "Multiple particles interaction test passed." << std::endl;
    } else {
        std::cerr << "Multiple particles interaction test failed." << std::endl;
    }
}

int main() {
    testDefaultConstructor();
    testParameterizedConstructor();
    testFullConstructor();
    testComputeDistance();
    testLennardJones();
    testCutLennardJones();
    testInteractionWithParticles();
    testMultipleParticles();
    return 0;
}
