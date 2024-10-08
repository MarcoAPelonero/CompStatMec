#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>

double adjustDelta(double (*distribution)(const std::vector<double>&), std::vector<double> start, int N, int M, double delta, double k);
void printResults(double result);
double getRandomValue(double min, double max);

#endif
