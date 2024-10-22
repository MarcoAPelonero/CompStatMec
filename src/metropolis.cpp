#include "metropolis.hpp"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>

// Constructor
Metropolis::Metropolis(double delta, int dimensions) : delta(delta), dimensions(dimensions) {}

// Multi-dimensional Metropolis integration
// Metropolis step in multiple dimensions
std::vector<double> Metropolis::metropolisStep(const std::vector<double>& x, double (*distribution)(const std::vector<double>&)) {
    std::vector<double> x_new = x;

    // Propose a new state by randomly perturbing each dimension
    for (int i = 0; i < dimensions; ++i) {
        x_new[i] = x[i] + (2.0 * delta * ((double)rand() / RAND_MAX - 0.5));

        // Apply boundary conditions: Ensure x_new stays within bounds [0, L]
        if (x_new[i] < 0.0) {
            x_new[i] = 0.0;
        } else if (x_new[i] > 100.0) {  // Assuming L = 100 for the map size
            x_new[i] = 100.0;
        }
    }

    // Compute the distribution values at the current and proposed points
    double current_value = distribution(x);
    double proposed_value = distribution(x_new);

    // Avoid division by zero - if current_value is zero, we reject unless proposed_value is non-zero
    if (current_value == 0.0) {
        if (proposed_value > 0.0) {
            return x_new;  // Accept if the new state is valid
        }
        return x;  // Reject if both are invalid (stay at current position)
    }

    // Compute the acceptance ratio
    double acceptanceRatio = proposed_value / current_value;

    // Accept or reject the new state based on the Metropolis criterion
    if (acceptanceRatio >= 0.9 || ((double)rand() / RAND_MAX) < acceptanceRatio) {
        return x_new;
    }

    return x;  // Reject the new state, stay at the current state
}


// Multi-dimensional Metropolis integration
double Metropolis::integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& start, int iterations) {
    std::vector<double> x = start;
    double sum = 0.0;
    int accepted_steps = 0;

    for (int i = 0; i < iterations; ++i) {
        std::vector<double> x_new = metropolisStep(x, distribution);
        if (x_new != x) {  // Count only accepted steps
            accepted_steps++;
        }
        sum += distribution(x_new);
        x = x_new;  // Move to the new position
    }
    
    return sum / accepted_steps;  // Normalize by accepted steps for a more accurate result
}

