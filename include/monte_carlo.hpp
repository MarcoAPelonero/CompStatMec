#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <vector>

class MonteCarlo {
public:
    MonteCarlo(int dimensions);
    
    // Classic Monte Carlo integration method
    double integrate(double (*distribution)(const std::vector<double>&), const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, int iterations);

private:
    int dimensions;  // Number of dimensions
};

#endif
