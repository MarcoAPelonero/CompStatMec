#include "monte_carlo_metropolis.hpp"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <numeric>

// This is the constructor method, and it's the __init__- for the cpp class
MonteCarloMetropolis::MonteCarloMetropolis(double delta, int dimensions) : delta(delta), dimensions(dimensions) {}

// This is the crucial method for the class, the integration part where basically the integral is evaluated and normalized before the return 
double MonteCarloMetropolis::integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& start, int iterations) {
    std::vector<double> x = start;
    double sum = 0.0;

    for (int i = 0; i < iterations; ++i) {
        std::vector<double> x_new = metropolisStep(x, distribution);
        sum += distribution(x_new);
        x = x_new;  
    }

    return sum / iterations;  
}

std::vector<double> MonteCarloMetropolis::metropolisStep(const std::vector<double>& x, double (*distribution)(const std::vector<double>&)) {
    std::vector<double> x_new = x;

    for (int i = 0; i < dimensions; ++i) {
        x_new[i] = x[i] + (2.0 * delta * ((double)rand() / RAND_MAX - 0.5));
    }

    double acceptanceRatio = distribution(x_new) / distribution(x);
    if (acceptanceRatio >= 1.0 || ((double)rand() / RAND_MAX) < acceptanceRatio) {
        return x_new;
    }

    return x;  
}

double MonteCarloMetropolis::classicMonteCarlo(double (*distribution)(const std::vector<double>&), const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, int iterations) {
    double sum = 0.0;
    std::vector<double> random_point(dimensions);

    for (int i = 0; i < iterations; ++i) {
        
        for (int j = 0; j < dimensions; ++j) {
            random_point[j] = lower_bound[j] + (upper_bound[j] - lower_bound[j]) * ((double)rand() / RAND_MAX);
        }
        sum += distribution(random_point);  
    }

    double volume = 1.0;
    for (int j = 0; j < dimensions; ++j) {
        volume *= (upper_bound[j] - lower_bound[j]);
    }

    return volume * sum / iterations;  
}

