#include "metropolis.hpp"
#include <cmath>
#include <cstdlib>
#include <vector>

// Constructor
Metropolis::Metropolis(double delta, int dimensions) : delta(delta), dimensions(dimensions) {}

// Multi-dimensional Metropolis integration
double Metropolis::integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& start, int iterations) {
    std::vector<double> x = start;
    double sum = 0.0;

    for (int i = 0; i < iterations; ++i) {
        std::vector<double> x_new = metropolisStep(x, distribution);
        sum += distribution(x_new);
        x = x_new;  // Move to the new position
    }

    return sum / iterations;
}

// Metropolis step in multiple dimensions
std::vector<double> Metropolis::metropolisStep(const std::vector<double>& x, double (*distribution)(const std::vector<double>&)) {
    std::vector<double> x_new = x;

    // Propose a new state by randomly perturbing each dimension
    for (int i = 0; i < dimensions; ++i) {
        x_new[i] = x[i] + (2.0 * delta * ((double)rand() / RAND_MAX - 0.5));
    }

    // Accept or reject the new state based on the Metropolis criterion
    double acceptanceRatio = distribution(x_new) / distribution(x);
    if (acceptanceRatio >= 1.0 || ((double)rand() / RAND_MAX) < acceptanceRatio) {
        return x_new;
    }

    return x;  // Reject the new state, stay at the current state
}