// FILE: CompStatMec/simulations/MCSim/src/test_particle.cpp

#include <iostream>
#include "particles.hpp"
#include "configuration.hpp"
#include "vec.hpp"

void testParticle() {
    std::cout << "Testing Particle class:" << std::endl;

    // Test default constructor
    Particle p1;
    std::cout << "Particle p1 (default constructor):" << std::endl;
    p1.getPosition().show("Position");

    // Test constructor with position vector
    Vector pos(1.0);
    Particle p2(pos);
    std::cout << "Particle p2 (position constructor):" << std::endl;
    p2.getPosition().show("Position");

    // Test constructor with position, velocity, acceleration, and mass
    Vector vel(2.0);
    Vector acc(3.0);
    ntype mass = 1.5;
    Particle p3(pos, vel, acc, mass);
    std::cout << "Particle p3 (full constructor):" << std::endl;
    p3.getPosition().show("Position");

    // Test store and restore
    p3.store();
    std::cout << "Particle p3 after random:" << std::endl;
    p3.getPosition().show("Position");
    p3.restore();
    std::cout << "Particle p3 after restore:" << std::endl;
    p3.getPosition().show("Position");
}

void testParticleEnsemble() {
    std::cout << "Testing particleEnsemble class:" << std::endl;

    // Test constructor with number of particles
    int N = 5;
    particleEnsemble ensemble(N, 10.0);
    std::cout << "particleEnsemble with " << N << " particles created." << std::endl;

    // Test initializeRandom
    ntype L = 10.0;
    ensemble.initializeRandom(L);
    std::cout << "particleEnsemble after initializeRandom:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Particle " << i + 1 << " position: ";
        ensemble(i).getPosition().show();
    }

    // Test store and restore
    ensemble.store();
    ensemble.initializeRandom(L);
    std::cout << "particleEnsemble after another initializeRandom:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Particle " << i + 1 << " position: ";
        ensemble(i).getPosition().show();
    }
    ensemble.restore();
    std::cout << "particleEnsemble after restore:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Particle " << i + 1 << " position: ";
        ensemble(i).getPosition().show();
    }
}

int main() {
    testParticle();
    testParticleEnsemble();
    return 0;
}
