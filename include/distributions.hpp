#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP
#include <vector>

double nile_distribution(const std::vector<double>& position);
double normalDistribution(double x);
double boltzmannDistribution(double x, double T);
#endif