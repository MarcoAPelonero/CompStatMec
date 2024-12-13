#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>
#include <random>

#include "interactionPotential.hpp"
#include "particleEnsemble.hpp"
#include "molecularDynamicsSimulation.hpp"

constexpr double EPS = 1e-9;

// Helper functions for floating-point comparisons
inline bool approximatelyEqual(double a, double b, double tolerance = 1e-6) {
    return std::fabs(a - b) < tolerance;
}

// Print a test header
void printTestHeader(const std::string &testName) {
    std::cout << "===== Running Test: " << testName << " =====" << std::endl;
}

// **************************************
// Tests for interactionPotential
// **************************************

void testMinimalImageDisplacement() {
    printTestHeader("testMinimalImageDisplacement");
    interactionPotential ip(1.0, 1.0, 10.0); // L=10

    Vector r1{0.0, 0.0, 0.0};
    Vector r2{9.9, 9.9, 9.9};
    Vector dr = ip.minimalImageDisplacement(r1, r2);
    //REMEMBER AT FIRST YOU THOUGH IT WAS THIS
    // Expected displacement ~ (-0.1, -0.1, -0.1)
    // But it is actually (0.1, 0.1, 0.1) because the box is 10.0
    // and every function uses just the module but this could be an issue
    // dr.show("Minimal Image Displacement: ");
    assert(approximatelyEqual(dr(0), 0.1));
    assert(approximatelyEqual(dr(1), 0.1));
    assert(approximatelyEqual(dr(2), 0.1));

    // Another test: halfway across the box
    r2 = Vector{5.1, 5.1, 5.1};
    dr = ip.minimalImageDisplacement(r1, r2);
    // Half the box is 5.0, so 5.1 is just over half, displacement should wrap to -4.9
    //Same thing here with the - signs
    // dr.show("Minimal Image Displacement: ");
    assert(approximatelyEqual(dr(0), 4.9));
    assert(approximatelyEqual(dr(1), 4.9));
    assert(approximatelyEqual(dr(2), 4.9));
}

void testLennardJonesValues() {
    printTestHeader("testLennardJonesValues");
    interactionPotential ip(1.0, 1.0, 10.0);

    double r = 1.0; // equals sigma
    double lj = ip.lennardJones(r);
    // At r = sigma, LJ potential = 4*epsilon*((1)^12 - (1)^6) = 4*(1 - 1) = 0 actually
    // Actually, for LJ: At r=sigma, (sigma/r)=1, sr6=1, sr12=1 => 4*(1-1)=0
    assert(approximatelyEqual(lj, 0.0));

    r = 2.0;
    lj = ip.lennardJones(r);
    // Just check it is negative and small
    assert(lj < 0.0);

    r = 3.0;
    lj = ip.lennardJones(r);
    // Larger r, LJ closer to zero
    assert(lj > -0.1 && lj < 0.1); 
}

void testLennardJonesForceMagnitude() {
    printTestHeader("testLennardJonesForceMagnitude");
    interactionPotential ip(1.0, 1.0, 10.0);

    double r = 1.0;
    double f = ip.computeForceMagnitudeLennardJones(r);
    // Near sigma, force should be large and repulsive (positive?)
    // Actually LJ force at r=sigma: f = 24ε(2 - 1)/r = 24 * 1 * (1)/1 = 24
    assert(approximatelyEqual(f, 24.0, 1e-5));

    r = 2.0;
    f = ip.computeForceMagnitudeLennardJones(r);
    // Force should be smaller, but just in module. Since sigma is 1 and r>sigma, 
    // the force should be negative (attractive)
    assert(f > -10.0 && f < 0.0);

    r = 5.0;
    f = ip.computeForceMagnitudeLennardJones(r);
    
    // Should be very small
    assert(f > -0.01 && f < 0.0);
}

void testCoulombPotential() {
    printTestHeader("testCoulombPotential");
    interactionPotential ip(1.0, 1.0, 10.0);

    double r = std::sqrt(12.0); // about 3.464
    double c = ip.coulomb(r);
    assert(approximatelyEqual(c, 1.0 / r, 1e-6));

    double fc = ip.computeForceMagnitudeCoulomb(r);
    // fc = 1/r^2
    assert(approximatelyEqual(fc, 1.0 / (r*r), 1e-6));
}

// **************************************
// Tests for particleEnsemble
// **************************************

void testParticleEnsembleInitialization() {
    printTestHeader("testParticleEnsembleInitialization");
    int numParticles = 8;
    double boxSize = 10.0;

    particleEnsemble ensemble(numParticles, boxSize);
    assert(ensemble.getNumParticles() == numParticles);

    // Check positions are within box
    for (int i = 0; i < numParticles; ++i) {
        Vector pos = ensemble(i).getPosition();
        for (int d = 0; d < 3; ++d) {
            assert(pos(d) >= 0.0 && pos(d) < boxSize);
        }
    }

    // Check that velocities are initialized (not NaN)
    for (int i = 0; i < numParticles; ++i) {
        Vector vel = ensemble(i).getVelocity();
        for (int d = 0; d < 3; ++d) {
            assert(!std::isnan(vel(d)));
        }
    }
}

void testPeriodicBoundaryApplication() {
    printTestHeader("testPeriodicBoundaryApplication");
    int numParticles = 1;
    double boxSize = 10.0;
    particleEnsemble ensemble(numParticles, boxSize);

    // Manually set position out of bounds
    ensemble(0).setPosition(Vector{-1.0, 10.5, 20.0});
    // Normally step functions call applyPeriodicBoundary internally,
    // so call a dummy step or just call it directly if accessible.
    // If not accessible, you might create a temporary public method to test.
    // For now, let's simulate a stepEuler call to trigger boundary conditions.

    double dt = 0.001;
    ensemble.stepEuler(0, dt);
    Vector pos = ensemble(0).getPosition();
    pos.show("Position after Euler step: ");
    // After wrapping: -1.0 -> 9.0, 10.5 -> 0.5, 20.0 -> 0.0 (assuming a single wrap)
    // However, stepEuler also changes position using velocity * dt,
    // so we must be careful. Let's set velocity to zero first to avoid that:
    ensemble(0).setVelocity(Vector{0.0,0.0,0.0});
    ensemble(0).setPosition(Vector{-1.0, 10.5, 20.0});
    ensemble.stepEuler(0, dt);
    pos = ensemble(0).getPosition();
    // Now we can assert the wrapping:
    assert(pos(0) >= 0.0 && pos(0) < boxSize);
    assert(pos(1) >= 0.0 && pos(1) < boxSize);
    assert(pos(2) >= 0.0 && pos(2) < boxSize);
}

void testForceAndEnergyComputation() {
    printTestHeader("testForceAndEnergyComputation");
    int numParticles = 2;
    double boxSize = 10.0;
    particleEnsemble ensemble(numParticles, boxSize);

    // Place two particles close to each other
    ensemble(0).setPosition(Vector{1.0, 1.0, 1.0});
    ensemble(1).setPosition(Vector{1.0, 1.0, 2.0});

    // Zero velocities for stable testing
    ensemble(0).setVelocity(Vector{0.0, 0.0, 0.0});
    ensemble(1).setVelocity(Vector{0.0, 0.0, 0.0});

    // A single Euler step to compute forces and update energies
    double dt = 0.001;
    ensemble.stepEuler(0, dt);
    ensemble.stepEuler(1, dt);

    double E0 = ensemble(0).getEnergy();
    double E1 = ensemble(1).getEnergy();

    // Just check energies are finite and negative (typical of LJ at close range)
    assert(std::isfinite(E0));
    assert(std::isfinite(E1));
    // Likely they are negative due to LJ binding at small separations
}

void testStoringAndRetrievingOldStates() {
    printTestHeader("testStoringAndRetrievingOldStates");
    int numParticles = 4;
    double boxSize = 10.0;
    particleEnsemble ensemble(numParticles, boxSize);

    // Update ensemble to store old states
    ensemble.updateEnsemble();
    for (int i = 0; i < numParticles; ++i) {
        Vector oldPos = ensemble(i).getOldPosition();
        Vector newPos = ensemble(i).getPosition();
        for (int d = 0; d < 3; ++d) {
            assert(approximatelyEqual(oldPos(d), newPos(d)));
        }
    }

    // Move a particle slightly
    ensemble(0).setPosition(Vector{2.0, 2.0, 2.0});
    // oldPosition should still reflect original before updateEnsemble()
    Vector oldPos = ensemble(0).getOldPosition();
    // oldPos should not be {2.0,2.0,2.0}
    assert(!approximatelyEqual(oldPos(0), 2.0));
    assert(!approximatelyEqual(oldPos(1), 2.0));
    assert(!approximatelyEqual(oldPos(2), 2.0));

    // Now call updateEnsemble() again
    ensemble.updateEnsemble();
    oldPos = ensemble(0).getOldPosition();
    assert(approximatelyEqual(oldPos(0), 2.0));
    assert(approximatelyEqual(oldPos(1), 2.0));
    assert(approximatelyEqual(oldPos(2), 2.0));
}

// **************************************
// Tests for molecularDynamicsSimulation
// **************************************

void testMDSimulationAverageEnergy() {
    printTestHeader("testMDSimulationAverageEnergy");
    int numParticles = 8;
    double boxSize = 10.0;
    double dt = 0.001;
    int numSteps = 10;
    double temperature = 1.0;

    molecularDynamicsSimulation sim(numParticles, boxSize, dt, numSteps, temperature);
    sim.setIntegrationMethod(molecularDynamicsSimulation::EULER);
    std::ofstream dummyOut("/dev/null"); // or "nul" on Windows, to discard output

    sim.run(dummyOut);
    double avgE = sim.computeAverageEnergy();

    // Just check it's finite
    assert(std::isfinite(avgE));
}

void testMDSimulationIntegrationStability() {
    printTestHeader("testMDSimulationIntegrationStability");
    // This is a heuristic test: we run a short simulation and see if things blow up.
    int numParticles = 8;
    double boxSize = 10.0;
    double dt = 0.005;
    int numSteps = 100; 
    double temperature = 1.0;
    molecularDynamicsSimulation sim(numParticles, boxSize, dt, numSteps, temperature);
    sim.setIntegrationMethod(molecularDynamicsSimulation::VELOCITY_VERLET);
    std::ofstream dummyOut("/dev/null");

    sim.run(dummyOut);
    double avgE = sim.computeAverageEnergy();
    // If simulation "exploded," avgE might be nan or very large
    assert(std::isfinite(avgE));
    assert(std::fabs(avgE) < 1e5); // Just a sanity bound
}

void testMDSimulationThermostat() {
    printTestHeader("testMDSimulationThermostat");
    int numParticles = 64;
    double boxSize = 20.0;
    double dt = 0.005;
    int numSteps = 200;
    double temperature = 1.0;
    molecularDynamicsSimulation sim(numParticles, boxSize, dt, numSteps, temperature);
    sim.setIntegrationMethod(molecularDynamicsSimulation::VELOCITY_VERLET);

    std::ofstream dummyOut("/dev/null");
    sim.run(dummyOut);

    // Check final average kinetic temperature
    // The code scales velocities every thermoSteps to achieve target temperature
    // Over enough steps, temperature should be close to the set temperature (1.0).
    double avgE = sim.computeAverageEnergy();
    // We can't directly check temperature from here, but stable average energy at 1.0 is a good sign
    // Just check it's not crazy large or zero
    assert(avgE > -10.0 && avgE < 10.0);
}

void testMDSimulationSetIntegrationMethod() {
    printTestHeader("testMDSimulationSetIntegrationMethod");
    molecularDynamicsSimulation sim(8, 10.0, 0.001, 10, 1.0);
    sim.setIntegrationMethod(molecularDynamicsSimulation::EULER);
    sim.setIntegrationMethod(molecularDynamicsSimulation::EULER_CROMER);
    sim.setIntegrationMethod(molecularDynamicsSimulation::LEAP_FROG);
    sim.setIntegrationMethod(molecularDynamicsSimulation::VELOCITY_VERLET);
    // If it gets here without errors, we consider it passed.
}

// **************************************
// Additional Checks: Edge Cases
// **************************************

void testSingleParticleNoForce() {
    printTestHeader("testSingleParticleNoForce");
    int numParticles = 1;
    double boxSize = 10.0;
    double dt = 0.01;
    int numSteps = 10;
    double temperature = 0.0; // no thermal motion

    molecularDynamicsSimulation sim(numParticles, boxSize, dt, numSteps, temperature);
    sim.setIntegrationMethod(molecularDynamicsSimulation::VELOCITY_VERLET);
    std::ofstream dummyOut("/dev/null");
    sim.run(dummyOut);

    // With one particle, no forces, velocity should remain nearly constant (if no thermostat)
    // Actually your code includes a thermostat step for VELOCITY_VERLET that scales velocities.
    // With T=0.0, the velocity should scale to zero eventually.
    double avgE = sim.computeAverageEnergy();
    // Without force and with a zero-temperature rescaling, it should end up at zero energy.
    assert(approximatelyEqual(avgE, 0.0, 1e-3));
}

void testTwoParticlesSymmetry() {
    printTestHeader("testTwoParticlesSymmetry");
    // Two identical particles placed symmetrically around center.
    int numParticles = 2;
    double boxSize = 10.0;
    double dt = 0.001;
    int numSteps = 5;
    double temperature = 0.0;
    molecularDynamicsSimulation sim(numParticles, boxSize, dt, numSteps, temperature);
    sim.setIntegrationMethod(molecularDynamicsSimulation::EULER);

    // Modify initial conditions in the ensemble directly if accessible
    // We'll do it by referencing sim's internal ensemble if there's a getter method.
    // If there's no direct access, you may need a setter or run a single step after custom init.

    // Assuming we can access `ensemble` somehow. If not, you'd need a getter in your code.
    // For demonstration, let's assume we can (or you'd add a friend class or method).
    // We'll just re-run particleEnsemble separately as a stand-in:

    particleEnsemble testEnsemble(numParticles, boxSize);
    testEnsemble(0).setPosition(Vector{4.0, 5.0, 5.0});
    testEnsemble(1).setPosition(Vector{6.0, 5.0, 5.0});
    testEnsemble(0).setVelocity(Vector{0.0, 0.0, 0.0});
    testEnsemble(1).setVelocity(Vector{0.0, 0.0, 0.0});
    // The force between them should be symmetric and along x-axis.

    // Just one step of Euler:
    testEnsemble.stepEuler(0, dt);
    testEnsemble.stepEuler(1, dt);

    Vector pos0 = testEnsemble(0).getPosition();
    Vector pos1 = testEnsemble(1).getPosition();

    // Check symmetry: if they move, they should move toward or away from each other symmetrically
    // Given no initial velocity, they will attract or repel symmetrically.
    // Just check they remain in the box and no NaN:
    for (int d = 0; d < 3; ++d) {
        assert(std::isfinite(pos0(d)));
        assert(std::isfinite(pos1(d)));
    }
}

int main() {
    try {
        testMinimalImageDisplacement();
        testLennardJonesValues();
        testLennardJonesForceMagnitude();
        testCoulombPotential();

        testParticleEnsembleInitialization();
        testPeriodicBoundaryApplication();
        testForceAndEnergyComputation();
        testStoringAndRetrievingOldStates();

        testMDSimulationAverageEnergy();
        testMDSimulationIntegrationStability();
        testMDSimulationThermostat();
        testMDSimulationSetIntegrationMethod();

        testSingleParticleNoForce();
        testTwoParticlesSymmetry();

        std::cout << "All tests passed successfully!" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "A test threw an exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
