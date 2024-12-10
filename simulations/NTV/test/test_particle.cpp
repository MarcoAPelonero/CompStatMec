#include <iostream>
#include "particle.hpp"

void testDefaultConstructor() {
    Particle p;
    std::cout << "Default constructor test passed." << std::endl;
}

void testParameterizedConstructor() {
    Vector r{1.032321, 2.021232, 3.04312321};
    Particle p(r);
    if (p.getPosition() == r) {
        std::cout << "Parameterized constructor test passed." << std::endl;
    } else {
        std::cerr << "Parameterized constructor test failed." << std::endl;
    }
}

void testFullConstructor() {
    Vector r{1.0, 2.0, 3.0};
    Vector v{0.121, 0.2123, 0.3524};
    Vector a{0.01, 0.02, 0.03};
    ntype m = 1.0;
    Particle p(r, v, a, m);
    if (p.getPosition() == r && p.getVelocity() == v && p.getAcceleration() == a) {
        std::cout << "Full constructor test passed." << std::endl;
    } else {
        std::cerr << "Full constructor test failed." << std::endl;
    }
}

void testStoreRestore() {
    Vector r{1.0, 2.0, 3.0};
    Particle p(r);
    p.store();
    Vector newR{4.0, 5.0, 6.0};
    p.setPosition(newR);
    p.restore();
    if (p.getPosition() == r) {
        std::cout << "Store and restore test passed." << std::endl;
    } else {
        std::cerr << "Store and restore test failed." << std::endl;
    }
}

void testRandom() {
    Particle p;
    p.random(10.0);
    Vector pos = p.getPosition();
    bool passed = true;
    for (int i = 0; i < dim; ++i) {
        if (pos.get(i) < 0.0 || pos.get(i) > 10.0) {
            passed = false;
            break;
        }
    }
    if (passed) {
        std::cout << "Random initialization test passed." << std::endl;
    } else {
        std::cerr << "Random initialization test failed." << std::endl;
    }
}

void testSetPosition() {
    Vector r{1.0, 2.0, 3.0};
    Particle p;
    p.setPosition(r);
    if (p.getPosition() == r) {
        std::cout << "Set position test passed." << std::endl;
    } else {
        std::cerr << "Set position test failed." << std::endl;
    }
}

void testSetVelocity() {
    Vector v{0.1, 0.2, 0.3};
    Particle p;
    p.setVelocity(v);
    if (p.getVelocity() == v) {
        std::cout << "Set velocity test passed." << std::endl;
    } else {
        std::cerr << "Set velocity test failed." << std::endl;
    }
}

void testAssignmentOperator() {
    Vector r{1.0, 2.0, 3.0};
    Particle p1(r);
    Particle p2;
    p2 = p1;
    if (p2.getPosition() == r) {
        std::cout << "Assignment operator test passed." << std::endl;
    } else {
        std::cerr << "Assignment operator test failed." << std::endl;
    }
}

int main() {
    testDefaultConstructor();
    testParameterizedConstructor();
    testFullConstructor();
    testStoreRestore();
    testRandom();
    testSetPosition();
    testSetVelocity();
    testAssignmentOperator();
    return 0;
}