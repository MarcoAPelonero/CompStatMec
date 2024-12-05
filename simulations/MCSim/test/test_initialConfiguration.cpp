// test_initialConfiguration.cpp

#include <iostream>
#include "initialConfiguration.hpp"
#include "configuration.hpp"

int main() {
    // Test default constructor
    initialConfiguration configDefault;
    particleEnsemble ensembleDefault = configDefault.returnParticleEnsemble();
    if (ensembleDefault.getNumParticles() == 0) {
        std::cout << "Default constructor test passed." << std::endl;
    } else {
        std::cout << "Default constructor test failed." << std::endl;
    }

    // Test parameterized constructor
    int N = 27;
    initialConfiguration configN(N);
    particleEnsemble ensembleN = configN.returnParticleEnsemble();
    if (ensembleN.getNumParticles() == N) {
        std::cout << "Parameterized constructor test passed." << std::endl;
    } else {
        std::cout << "Parameterized constructor test failed." << std::endl;
    }

    // Test initializeSquareLattice
    ntype L = 10.0;
    configN.initializeSquareLattice(L);
    ensembleN = configN.returnParticleEnsemble();
    bool latticeTestPassed = true;

    // Check that particles are within the box of size L
    for (int i = 0; i < N; ++i) {
        Vector pos = ensembleN(i).getPosition();
        for (int cord; cord < 3; cord++) {
            if (pos(cord) < 0 || pos(cord) > L) {
                latticeTestPassed = false;
                break;
            }
        }
    }

    if (latticeTestPassed) {
        std::cout << "initializeSquareLattice test passed." << std::endl;
    } else {
        std::cout << "initializeSquareLattice test failed." << std::endl;
    }

    // Additional tests for initializeRandom can be added similarly

    return 0;
}