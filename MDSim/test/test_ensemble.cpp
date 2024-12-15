#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include "vec.hpp"
#include "particle.hpp"
#include "particleEnsemble.hpp"

bool areVectorsEqual(const Vec &v1, const Vec &v2, double tolerance = 1e-6) {
    return (std::fabs(v1.getX() - v2.getX()) < tolerance &&
            std::fabs(v1.getY() - v2.getY()) < tolerance &&
            std::fabs(v1.getZ() - v2.getZ()) < tolerance);
}

void testInitialization() {
    int numParticles = 27;
    double boxLength = 3.0;
    ParticleEnsemble ensemble(numParticles, boxLength);

    const std::vector<Particle>& particles = ensemble.getParticles();
    bool noOverlap = true;
    double expectedSpacing = boxLength / std::ceil(std::cbrt(numParticles));

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Vec pos1 = particles[i].getPosition();
            Vec pos2 = particles[j].getPosition();
            double distance = std::sqrt(std::pow(pos1.getX() - pos2.getX(), 2) +
                                        std::pow(pos1.getY() - pos2.getY(), 2) +
                                        std::pow(pos1.getZ() - pos2.getZ(), 2));
            if (distance < expectedSpacing) {
                noOverlap = false;
                break;
            }
        }
    }

    if (noOverlap) {
        std::cout << "Initialization test passed: No overlap detected." << std::endl;
    } else {
        std::cout << "Initialization test failed: Overlap detected." << std::endl;
    }
}

void testVelocityDistribution() {
    int numParticles = 10000;
    double boxLength = 10.0;
    ParticleEnsemble ensemble(numParticles, boxLength);

    const std::vector<Particle>& particles = ensemble.getParticles();
    double sumX = 0, sumY = 0, sumZ = 0;
    double sumX2 = 0, sumY2 = 0, sumZ2 = 0;

    for (const Particle& p : particles) {
        Vec v = p.getVelocity();
        sumX += v.getX();
        sumY += v.getY();
        sumZ += v.getZ();
        sumX2 += v.getX() * v.getX();
        sumY2 += v.getY() * v.getY();
        sumZ2 += v.getZ() * v.getZ();
    }

    double meanX = sumX / numParticles;
    double meanY = sumY / numParticles;
    double meanZ = sumZ / numParticles;
    double stdX = std::sqrt(sumX2 / numParticles - meanX * meanX);
    double stdY = std::sqrt(sumY2 / numParticles - meanY * meanY);
    double stdZ = std::sqrt(sumZ2 / numParticles - meanZ * meanZ);

    std::cout << "Velocity distribution test results:" << std::endl;
    std::cout << "Mean X: " << meanX << ", Mean Y: " << meanY << ", Mean Z: " << meanZ << std::endl;
    std::cout << "Std X: " << stdX << ", Std Y: " << stdY << ", Std Z: " << stdZ << std::endl;

    if (std::abs(meanX) < 0.1 && std::abs(meanY) < 0.1 && std::abs(meanZ) < 0.1 &&
        std::abs(stdX - 1.0) < 0.1 && std::abs(stdY - 1.0) < 0.1 && std::abs(stdZ - 1.0) < 0.1) {
        std::cout << "Velocity distribution test passed." << std::endl;
    } else {
        std::cout << "Velocity distribution test failed." << std::endl;
    }
}
// Function to compute the analytical Lennard-Jones force
Vec analyticalLennardJonesForce(const Vec& r, double boxLength) {
    // Apply minimum image convention
    double x = r.getX();
    double y = r.getY();
    double z = r.getZ();

    if (x > boxLength / 2) x -= boxLength;
    if (x <= -boxLength / 2) x += boxLength;
    if (y > boxLength / 2) y -= boxLength;
    if (y <= -boxLength / 2) y += boxLength;
    if (z > boxLength / 2) z -= boxLength;
    if (z <= -boxLength / 2) z += boxLength;

    Vec r_min(x, y, z);
    double r_len = r_min.length();

    if (r_len == 0.0 || r_len > 4.0) {
        return Vec(0.0, 0.0, 0.0);
    }

    double r2 = r_len * r_len;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    double force_magnitude = 24.0 * (2.0 / r12 - 1.0 / r6) / r2;

    return Vec(force_magnitude * r_min.getX(),
               force_magnitude * r_min.getY(),
               force_magnitude * r_min.getZ());
}

void testLennardJonesForce() {
    double boxLength = 10.0; // Arbitrary box length
    ParticleEnsemble ensemble(0, boxLength); // Assuming constructor takes numParticles and boxLength

    // Ensure the ensemble starts empty
    ensemble.getParticles().clear();

    // Define test cases with pairs of particles
    struct TestCase {
        Vec position1;
        Vec position2;
        std::string description;
    };

    // Define multiple test cases
    TestCase testCases[] = {
        { Vec(0.0, 0.0, 0.0), Vec(0.5, 0.0, 0.0), "Strong Repulsion along +x" },
        { Vec(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), "Moderate Repulsion along +x" },
        { Vec(0.0, 0.0, 0.0), Vec(1.12246, 0.0, 0.0), "Equilibrium Distance along +x" },
        { Vec(0.0, 0.0, 0.0), Vec(2.0, 0.0, 0.0), "Mild Attraction along +x" },
        { Vec(0.0, 0.0, 0.0), Vec(4.5, 0.0, 0.0), "Beyond Cutoff along +x" },
        { Vec(1.0, 1.0, 1.0), Vec(9.0, 1.0, 1.0), "Minimum Image across -x boundary" },
        { Vec(5.0, 5.0, 5.0), Vec(5.0, 5.0, 5.0), "Identical Positions (Overlap)" },
        { Vec(2.0, 3.0, 4.0), Vec(7.0, 3.0, 4.0), "Minimum Image across -x boundary at distance 5.0" },
        { Vec(1.0, 1.0, 1.0), Vec(2.0, 2.0, 2.0), "Force at arbitrary angle" },
    };

    for (const auto& test : testCases) {
        // Create particles
        Particle p1(1.0, test.position1, Vec(0.0, 0.0, 0.0));
        Particle p2(1.0, test.position2, Vec(0.0, 0.0, 0.0));

        // Add particles to the ensemble
        std::vector<Particle>& particles = ensemble.getParticles();
        particles.clear();
        particles.push_back(p1);
        particles.push_back(p2);

        // Compute forces
        ensemble.computeForces();

        // Extract forces
        Vec forceOnP1 = particles[0].getForce();
        Vec forceOnP2 = particles[1].getForce();

        // Compute the displacement vector
        Vec displacement = particles[0].getPosition() - particles[1].getPosition();

        // Compute analytical force
        Vec expectedForce = analyticalLennardJonesForce(displacement, boxLength);

        // Print test description
        std::cout << "Test Case: " << test.description << "\n";
        std::cout << "Positions:\n";
        std::cout << "  Particle 1: (" << particles[0].getPosition().getX() << ", " 
                  << particles[0].getPosition().getY() << ", " << particles[0].getPosition().getZ() << ")\n";
        std::cout << "  Particle 2: (" << particles[1].getPosition().getX() << ", " 
                  << particles[1].getPosition().getY() << ", " << particles[1].getPosition().getZ() << ")\n";

        // Print computed forces
        std::cout << "Computed Forces:\n";
        std::cout << "  Force on Particle 1: (" << forceOnP1.getX() << ", " 
                  << forceOnP1.getY() << ", " << forceOnP1.getZ() << ")\n";
        std::cout << "  Force on Particle 2: (" << forceOnP2.getX() << ", " 
                  << forceOnP2.getY() << ", " << forceOnP2.getZ() << ")\n";

        // Print expected forces
        std::cout << "Expected Force on Particle 1: (" << expectedForce.getX() << ", " 
                  << expectedForce.getY() << ", " << expectedForce.getZ() << ")\n";
        std::cout << "Expected Force on Particle 2: (" << -expectedForce.getX() << ", " 
                  << -expectedForce.getY() << ", " << -expectedForce.getZ() << ")\n";

        // Calculate differences
        double diffFx = std::abs(forceOnP1.getX() - expectedForce.getX());
        double diffFy = std::abs(forceOnP1.getY() - expectedForce.getY());
        double diffFz = std::abs(forceOnP1.getZ() - expectedForce.getZ());

        double tolerance = 1e-5; // Adjust as needed

        // Validate the force on Particle 1
        assert(diffFx < tolerance && "Force X-component on Particle 1 is incorrect.");
        assert(diffFy < tolerance && "Force Y-component on Particle 1 is incorrect.");
        assert(diffFz < tolerance && "Force Z-component on Particle 1 is incorrect.");

        // Validate the force on Particle 2 (should be opposite)
        assert(std::abs(forceOnP2.getX() + expectedForce.getX()) < tolerance && 
               "Force X-component on Particle 2 is not opposite to Particle 1.");
        assert(std::abs(forceOnP2.getY() + expectedForce.getY()) < tolerance && 
               "Force Y-component on Particle 2 is not opposite to Particle 1.");
        assert(std::abs(forceOnP2.getZ() + expectedForce.getZ()) < tolerance && 
               "Force Z-component on Particle 2 is not opposite to Particle 1.");

        // Special checks based on description
        if (test.description.find("Overlap") != std::string::npos) {
            // At overlapping positions, force should be zero to avoid infinite forces
            assert(forceOnP1.length() == 0.0 && "Force should be zero for overlapping particles.");
            assert(forceOnP2.length() == 0.0 && "Force should be zero for overlapping particles.");
            std::cout << "  Overlapping particles: Forces correctly set to zero.\n";
        } else if (test.description.find("Minimum Image") != std::string::npos) {
            // Ensure that minimum image convention was applied correctly
            // e.g., for particles across the boundary, displacement should be adjusted
            // Already handled in analytical force
            // Just ensure that forces are computed correctly
            // Already covered by analytical force comparison
        }

        // Check Newton's Third Law: F12 = -F21
        double totalFx = forceOnP1.getX() + forceOnP2.getX();
        double totalFy = forceOnP1.getY() + forceOnP2.getY();
        double totalFz = forceOnP1.getZ() + forceOnP2.getZ();

        assert(std::abs(totalFx) < tolerance && "Newton's Third Law violated in X-component.");
        assert(std::abs(totalFy) < tolerance && "Newton's Third Law violated in Y-component.");
        assert(std::abs(totalFz) < tolerance && "Newton's Third Law violated in Z-component.");

        std::cout << "  Forces validated successfully.\n";
        std::cout << "----------------------------------------\n";
    }

    std::cout << "All Lennard-Jones force tests passed successfully.\n";
}

int main() {
    testInitialization();
    testVelocityDistribution();
    testLennardJonesForce();
    return 0;
}