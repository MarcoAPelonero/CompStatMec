#include "distributions.hpp"
#include <cmath>

double normalDistribution(double x) {
    return exp(-x * x / 2.0) / sqrt(2 * M_PI);
}

double boltzmannDistribution(double x, double T) {
    return exp(-x / T);
}
