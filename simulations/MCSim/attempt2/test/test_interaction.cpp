#include <iostream>
#include <cmath>
#include <cassert>

// Include your library headers
#include "interactionPotential.hpp"
#include "vec.hpp"
#include "particle.hpp"

// Define a small tolerance for floating-point comparisons
const double TOL = 1e-6;

// Helper function for floating-point comparison
bool almostEqual(double a, double b, double tol = TOL) {
    return std::abs(a - b) < tol;
}

// ==================== interactionPotential Class Tests ====================

void testInteractionPotentialConstructors() {
    // Default constructor
    interactionPotential ip1;
    assert(almostEqual(ip1.getSigma(), 1.0));
    assert(almostEqual(ip1.getEpsilon(), 1.0));
    assert(almostEqual(ip1.getL(), 10.0));
    assert(almostEqual(ip1.getRcut(), 2.5));
    std::cout << "interactionPotential default constructor test passed.\n";

    // Parameterized constructor (box length)
    interactionPotential ip2(15.0);
    assert(almostEqual(ip2.getL(), 15.0));
    assert(almostEqual(ip2.getSigma(), 1.0));
    assert(almostEqual(ip2.getEpsilon(), 1.0));
    assert(almostEqual(ip2.getRcut(), 2.5));
    std::cout << "interactionPotential parameterized constructor (box length) test passed.\n";

    // Full parameterized constructor
    interactionPotential ip3(1.5, 0.8, 20.0, 3.0);
    assert(almostEqual(ip3.getSigma(), 1.5));
    assert(almostEqual(ip3.getEpsilon(), 0.8));
    assert(almostEqual(ip3.getL(), 20.0));
    assert(almostEqual(ip3.getRcut(), 3.0));
    std::cout << "interactionPotential full parameterized constructor test passed.\n";
}

void testComputeDistance() {
    interactionPotential ip;
    Vector r1{1.0, 1.0, 1.0};
    Vector r2{4.0, 5.0, 6.0};
    double distance = ip.computeDistance(r1, r2);
    double expected = std::sqrt( (1.0-4.0)*(1.0-4.0) + 
                                 (1.0-5.0)*(1.0-5.0) + 
                                 (1.0-6.0)*(1.0-6.0) );
    assert(almostEqual(distance, expected));
    std::cout << "interactionPotential computeDistance test without MIC passed.\n";

    // Test with MIC
    // Assume box length L=10.0, r1=(1,1,1), r2=(9,9,9)
    Vector r3{9.0, 9.0, 9.0};
    double distance_mic = ip.computeDistance(r1, r3);
    // Minimum image: dr = (1-9) = -8 => since |dr| > L/2=5, dr += L => dr = -8 + 10 = 2
    // Thus distance = sqrt(3*(2^2)) = sqrt(12) ≈ 3.4641
    double expected_mic = std::sqrt(12.0);
    assert(almostEqual(distance_mic, expected_mic));
    std::cout << "interactionPotential computeDistance test with MIC passed.\n";
}

void testLennardJones() {
    interactionPotential ip;
    double r = 1.0;
    double lj = ip.lennardJones(r);
    double expectedLJ = 4.0 * ip.getEpsilon() * (std::pow(ip.getSigma() / r, 12) - std::pow(ip.getSigma() / r, 6));
    assert(almostEqual(lj, expectedLJ));
    std::cout << "interactionPotential Lennard-Jones potential test passed.\n";

    // Test at different distance
    r = 2.0;
    lj = ip.lennardJones(r);
    expectedLJ = 4.0 * ip.getEpsilon() * (std::pow(ip.getSigma() / r, 12) - std::pow(ip.getSigma() / r, 6));
    assert(almostEqual(lj, expectedLJ));
    std::cout << "interactionPotential Lennard-Jones potential at r=2.0 test passed.\n";
}

void testCutLennardJones() {
    interactionPotential ip;
    double r_inside = 2.0;
    double r_outside = 3.0;

    double clj_inside = ip.cutLennardJones(r_inside);
    double expected_clj_inside = ip.lennardJones(r_inside) - ip.lennardJones(ip.getRcut());
    assert(almostEqual(clj_inside, expected_clj_inside));
    std::cout << "interactionPotential cut Lennard-Jones potential inside cutoff test passed.\n";

    // Exactly at cutoff
    double clj_cutoff = ip.cutLennardJones(ip.getRcut());
    // Depending on implementation, at r = rcut, it might be zero or V(rcut) - V(rcut) = 0
    assert(almostEqual(clj_cutoff, 0.0));
    std::cout << "interactionPotential cut Lennard-Jones potential at cutoff distance test passed.\n";

    // Outside cutoff
    double clj_outside = ip.cutLennardJones(r_outside);
    double expected_clj_outside = 0.0;
    assert(almostEqual(clj_outside, expected_clj_outside));
    std::cout << "interactionPotential cut Lennard-Jones potential outside cutoff test passed.\n";
}

void testCutLennardJonesParameters() {
    // Test with different sigma and epsilon
    interactionPotential ip(2.0, 0.5, 10.0, 3.0);
    double r = 2.0;
    double lj = ip.lennardJones(r);
    double expectedLJ = 4.0 * ip.getEpsilon() * (std::pow(ip.getSigma() / r, 12) - std::pow(ip.getSigma() / r, 6));
    assert(almostEqual(lj, expectedLJ));
    std::cout << "interactionPotential Lennard-Jones potential with custom parameters test passed.\n";

    double clj = ip.cutLennardJones(r);
    double expected_clj = lj - ip.lennardJones(ip.getRcut());
    assert(almostEqual(clj, expected_clj));
    std::cout << "interactionPotential cut Lennard-Jones potential with custom parameters inside cutoff test passed.\n";
}

void testComputeDistancePeriodicBoundaries() {
    interactionPotential ip(1.0, 1.0, 10.0, 2.5); // Using default sigma and epsilon
    Vector r1{1.0, 1.0, 1.0};
    Vector r2{9.0, 9.0, 9.0}; // Should use MIC

    double distance = ip.computeDistance(r1, r2);
    double expected_distance = std::sqrt(3 * std::pow(2.0, 2)); // dr = -8 + 10 = 2 for each component
    expected_distance = std::sqrt(12.0); // ≈ 3.4641
    assert(almostEqual(distance, expected_distance));
    std::cout << "interactionPotential computeDistance with periodic boundaries (MIC) test passed.\n";
}

void testMultipleInteractions() {
    interactionPotential ip;
    // Define multiple particle positions
    std::vector<Vector> positions = {
        Vector{1.0, 1.0, 1.0},
        Vector{2.0, 2.0, 2.0},
        Vector{3.0, 3.0, 3.0},
        Vector{8.0, 8.0, 8.0} // Should interact with particle 0 via MIC
    };

    // Expected distances
    std::vector<std::vector<double>> expected_distances = {
        {0.0, std::sqrt(3.0), std::sqrt(12.0), std::sqrt(12.0)}, // Particle 0
        {std::sqrt(3.0), 0.0, std::sqrt(3.0), std::sqrt(12.0)},   // Particle 1
        {std::sqrt(12.0), std::sqrt(3.0), 0.0, std::sqrt(27.0)}, // Particle 2
        {std::sqrt(12.0), std::sqrt(12.0), std::sqrt(27.0), 0.0} // Particle 3
    };

    // Compute and verify distances
    for(std::vector<Vector>::size_type i = 0; i < positions.size(); ++i) {
        for(std::vector<Vector>::size_type j = 0; j < positions.size(); ++j) {
            double computed_distance = ip.computeDistance(positions[i], positions[j]);
            if(i == j) {
                assert(almostEqual(computed_distance, 0.0));
            } else {
                // Apply MIC manually
                Vector dr = positions[i] - positions[j];
                for(int d = 0; d < dim; ++d) {
                    if(dr.get(d) > ip.getL() / 2.0) {
                        dr.set(d, dr.get(d) - ip.getL());
                    }
                    else if(dr.get(d) < -ip.getL() / 2.0) {
                        dr.set(d, dr.get(d) + ip.getL());
                    }
                }
                double expected_distance_mic = dr.modulus();
                assert(almostEqual(computed_distance, expected_distance_mic));
            }
        }
    }
    std::cout << "interactionPotential multiple interactions distance computations test passed.\n";
}

void testEnergyCalculationConsistency() {
    // Create two particles and verify energy calculations
    interactionPotential ip;
    Vector r1{1.0, 1.0, 1.0};
    Vector r2{2.0, 2.0, 2.0};

    double distance = ip.computeDistance(r1, r2);
    double lj = ip.lennardJones(distance);
    double clj = ip.cutLennardJones(distance);

    // Manual calculation
    double expected_lj = 4.0 * ip.getEpsilon() * (std::pow(ip.getSigma() / distance, 12) - std::pow(ip.getSigma() / distance, 6));
    double expected_clj = (distance < ip.getRcut()) ? (expected_lj - 4.0 * ip.getEpsilon() * (std::pow(ip.getSigma() / ip.getRcut(), 12) - std::pow(ip.getSigma() / ip.getRcut(), 6))) : 0.0;

    assert(almostEqual(lj, expected_lj));
    assert(almostEqual(clj, expected_clj));

    std::cout << "interactionPotential energy calculation consistency test passed.\n";
}

void testEdgeCases() {
    interactionPotential ip;
    // Zero distance (should handle gracefully, possibly by returning a very large potential)
    Vector r1{1.0, 1.0, 1.000000000000001};
    Vector r2{1.0, 1.0, 1.0};
    double distance = ip.computeDistance(r1, r2);
    assert(almostEqual(distance, 0.0));
    double lj = ip.lennardJones(distance);
    // Depending on implementation, this might result in infinity or a very large number
    // Here, we check if it's greater than a threshold
    assert(lj > 1e6);
    std::cout << "interactionPotential zero distance edge case test passed.\n";

    double clj = ip.cutLennardJones(ip.getRcut());
    assert(almostEqual(clj, 0.0));
    std::cout << "interactionPotential distance exactly at cutoff test passed.\n";
}

// ==================== Main Test Runner ====================

int main() {
    std::cout << "Starting interactionPotential class tests...\n\n";

    testInteractionPotentialConstructors();
    testComputeDistance();
    testComputeDistancePeriodicBoundaries();
    testLennardJones();
    testCutLennardJones();
    testCutLennardJonesParameters();
    testMultipleInteractions();
    testEnergyCalculationConsistency();
    testEdgeCases();

    std::cout << "\nAll interactionPotential tests passed successfully!\n";
    return 0;
}

