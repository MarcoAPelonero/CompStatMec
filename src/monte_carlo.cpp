#include "monte_carlo.hpp"
#include <cstdlib>
#include <vector>

// Constructor
MonteCarlo::MonteCarlo(int dimensions) : dimensions(dimensions) {}

// Classic Monte Carlo integration over multiple dimensions
double MonteCarlo::integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, int iterations) {
    double sum = 0.0;
    std::vector<double> random_point(dimensions);

    for (int i = 0; i < iterations; ++i) {
        // Generate a random point in the defined multi-dimensional space
        for (int j = 0; j < dimensions; ++j) {
            random_point[j] = lower_bound[j] + (upper_bound[j] - lower_bound[j]) * ((double)rand() / RAND_MAX);
        }
        sum += distribution(random_point);  // Evaluate the function at the random point
    }

    // Compute the volume of the integration region
    double volume = 1.0;
    for (int j = 0; j < dimensions; ++j) {
        volume *= (upper_bound[j] - lower_bound[j]);
    }

    return volume * sum / iterations;  // Scale by the volume of the region
}
