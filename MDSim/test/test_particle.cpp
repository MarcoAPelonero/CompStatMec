#include <iostream>
#include <cmath>
#include "vec.hpp"
#include "particle.hpp"

void testDefaultConstructor() {
    Particle p;
    if (p.getMass() == 0 && p.getPosition() == Vec() && p.getVelocity() == Vec() && p.getForce() == Vec() &&
        p.getKineticEnergy() == 0 && p.getPotentialEnergy() == 0 && p.getTotalEnergy() == 0) {
        std::cout << "Default constructor test passed." << std::endl;
    } else {
        std::cout << "Default constructor test failed." << std::endl;
    }
}

void testParameterizedConstructor() {
    Vec position(1.0, 2.0, 3.0);
    Vec velocity(4.0, 5.0, 6.0);
    Particle p(10.0, position, velocity);
    if (p.getMass() == 10.0 && p.getPosition() == position && p.getVelocity() == velocity && p.getForce() == Vec() &&
        p.getKineticEnergy() == 0 && p.getPotentialEnergy() == 0 && p.getTotalEnergy() == 0) {
        std::cout << "Parameterized constructor test passed." << std::endl;
    } else {
        std::cout << "Parameterized constructor test failed." << std::endl;
    }
}

void testSettersAndGetters() {
    Particle p;
    Vec position(1.0, 2.0, 3.0);
    Vec velocity(4.0, 5.0, 6.0);
    Vec force(7.0, 8.0, 9.0);
    p.setMass(10.0);
    p.setPosition(position);
    p.setVelocity(velocity);
    p.setForce(force);
    p.setKineticEnergy(20.0);
    p.setPotentialEnergy(30.0);
    p.setTotalEnergy(50.0);

    if (p.getMass() == 10.0 && p.getPosition() == position && p.getVelocity() == velocity && p.getForce() == force &&
        p.getKineticEnergy() == 20.0 && p.getPotentialEnergy() == 30.0 && p.getTotalEnergy() == 50.0) {
        std::cout << "Setters and getters test passed." << std::endl;
    } else {
        std::cout << "Setters and getters test failed." << std::endl;
    }
}

void testComputeKineticEnergy() {
    Vec velocity(3.0, 4.0, 0.0);
    Particle p(2.0, Vec(), velocity);
    p.computeKineticEnergy();
    if (std::abs(p.getKineticEnergy() - 25.0) < 1e-9) {
        std::cout << "Compute kinetic energy test passed." << std::endl;
    } else {
        std::cout << "Compute kinetic energy test failed." << std::endl;
    }
}

void testComputeTotalEnergy() {
    Particle p;
    p.setKineticEnergy(20.0);
    p.setPotentialEnergy(30.0);
    p.computeTotalEnergy();
    if (std::abs(p.getTotalEnergy() - 50.0) < 1e-9) {
        std::cout << "Compute total energy test passed." << std::endl;
    } else {
        std::cout << "Compute total energy test failed." << std::endl;
    }
}

void testShow() {
    Vec position(1.0, 2.0, 3.0);
    Vec velocity(4.0, 5.0, 6.0);
    Vec force(7.0, 8.0, 9.0);
    Particle p(10.0, position, velocity);
    p.setForce(force);
    p.setKineticEnergy(20.0);
    p.setPotentialEnergy(30.0);
    p.setTotalEnergy(50.0);

    std::cout << "Expected output:" << std::endl;
    std::cout << "Mass: 10" << std::endl;
    std::cout << "Position: (1, 2, 3)" << std::endl;
    std::cout << "Velocity: (4, 5, 6)" << std::endl;
    std::cout << "Force: (7, 8, 9)" << std::endl;
    std::cout << "Kinetic energy: 20" << std::endl;
    std::cout << "Potential energy: 30" << std::endl;
    std::cout << "Total energy: 50" << std::endl;

    std::cout << "Actual output:" << std::endl;
    p.show();
}

int main() {
    testDefaultConstructor();
    testParameterizedConstructor();
    testSettersAndGetters();
    testComputeKineticEnergy();
    testComputeTotalEnergy();
    testShow();
    return 0;
}