#ifndef MONTE_CARLO_METROPOLIS_HPP
#define MONTE_CARLO_METROPOLIS_HPP

#include <vector>

class MonteCarloMetropolis {
public:
    MonteCarloMetropolis(double delta, int dimensions);
    
    // Multi-dimensional integration using the Metropolis algorithm
    double integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& start, int iterations);

    // Classic Monte Carlo integration method
    double classicMonteCarlo(double (*distribution)(const std::vector<double>&), const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, int iterations);

private:
    double delta;  // Step size for the Metropolis method
    int dimensions;  // Number of dimensions
    std::vector<double> metropolisStep(const std::vector<double>& x, double (*distribution)(const std::vector<double>&));
};

#endif

