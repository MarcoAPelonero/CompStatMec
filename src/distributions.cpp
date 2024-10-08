#include "distributions.hpp"
#include <cmath>

double nile_distribution(const std::vector<double>& position) {
    int L = 100;  // Map size (LxL)
    int nile_width = 10;  // Nile's width
    int nile_height = 20;  // Nile's height
    
    double x = position[0];
    double y = position[1];
    
    if (x >= (L/2 - nile_width/2) && x <= (L/2 + nile_width/2) &&
        y >= (L/2 - nile_height/2) && y <= (L/2 + nile_height/2)) {
        return 1.0;  // Inside the Nile
    }
    return 0.0;  // Outside the Nile
}

double normalDistribution(double x) {
    return exp(-x * x / 2.0) / sqrt(2 * M_PI);
}

double boltzmannDistribution(double x, double T) {
    return exp(-x / T);
}