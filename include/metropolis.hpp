#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include <vector>

class Metropolis {
public:
    Metropolis(double delta, int dimensions);

    // Multi-dimensional integration using the Metropolis algorithm
    double integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& start, int iterations);

private:
    double delta;  // Step size
    int dimensions;  // Number of dimensions
    std::vector<double> metropolisStep(const std::vector<double>& x, double (*distribution)(const std::vector<double>&));
};

#endif