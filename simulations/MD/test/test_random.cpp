#include <iostream>
#include <vector>
#include "randNumGen.hpp"

void testSeed() {
    randNumGen<> rng1, rng2;
    rng1.seed(42);
    rng2.seed(42);
    bool passed = true;
    for (int i = 0; i < 10; ++i) {
        if (rng1() != rng2()) {
            passed = false;
            break;
        }
    }
    if (passed) {
        std::cout << "Seed test passed." << std::endl;
    } else {
        std::cerr << "Seed test failed." << std::endl;
    }
}

void testRSeed() {
    randNumGen<> rng1, rng2;
    rng1.rseed();
    rng2.rseed();
    bool passed = true;
    for (int i = 0; i < 10; ++i) {
        if (rng1() == rng2()) {
            passed = false;
            break;
        }
    }
    if (passed) {
        std::cout << "RSeed test passed." << std::endl;
    } else {
        std::cerr << "RSeed test failed." << std::endl;
    }
}

void testSetDistribution() {
    randNumGen<> rng;
    rng.setDistribution(10.0, 20.0);
    bool passed = true;
    for (int i = 0; i < 10; ++i) {
        double num = rng();
        if (num < 10.0 || num > 20.0) {
            passed = false;
            break;
        }
    }
    if (passed) {
        std::cout << "Set distribution test passed." << std::endl;
    } else {
        std::cerr << "Set distribution test failed." << std::endl;
    }
}

void testMultipleGenerators() {
    const int numGenerators = 23;
    std::vector<randNumGen<>> generators1(numGenerators);
    std::vector<randNumGen<>> generators2(numGenerators);
    for (int i = 0; i < numGenerators; ++i) {
        generators1[i].seed(42);
        generators2[i].seed(42);
    }
    bool passed = true;
    for (int i = 0; i < numGenerators; ++i) {
        for (int j = 0; j < 10; ++j) {
            if (generators1[i]() != generators2[i]()) {
                passed = false;
                break;
            }
        }
        if (!passed) break;
    }
    if (passed) {
        std::cout << "Multiple generators with same seed test passed." << std::endl;
    } else {
        std::cerr << "Multiple generators with same seed test failed." << std::endl;
    }
}

int main() {
    testSeed();
    testRSeed();
    testSetDistribution();
    testMultipleGenerators();
    return 0;
}