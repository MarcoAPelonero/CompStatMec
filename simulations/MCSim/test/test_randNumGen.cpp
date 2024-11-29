#include <iostream>
#include "../include/randNumGen.hpp"

int main() {
    // Test 1: Initialize and seed the random number generator
    std::cout << "Testing randNumGen..." << std::endl;

    randNumGen<double, std::mt19937_64> rng;

    std::cout << "Seeding RNG with 42..." << std::endl;
    rng.seed(42);

    // Test 2: Generate some random numbers
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::cout << "Random numbers in [0, 1]: " << std::endl;
    for (int i = 0; i < 5; ++i) {
        double num = dist(rng);
        std::cout << "Random number " << i + 1 << ": " << num << std::endl;
    }

    // Success message
    std::cout << "Tests completed successfully!" << std::endl;

    return 0;
}