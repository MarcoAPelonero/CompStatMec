#include "utils.hpp"
#include <iostream>
#include <cstdlib>
#include <cmath>

void printResults(double result) {
    std::cout << "Integration result: " << result << std::endl;
}

double getRandomValue(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

double adjustDelta(double (*distribution)(const std::vector<double>&), std::vector<double> start, int N, int M, double delta, double k) {
    int accepted = 0;
    std::vector<double> x = start;
    int dimensions = start.size();
    
    // Simulate N steps
    for (int i = 0; i < N; ++i) {
        std::vector<double> x_new = x;

        // Propose new state in each dimension
        for (int j = 0; j < dimensions; ++j) {
            x_new[j] = x[j] + (2.0 * delta * ((double)rand() / RAND_MAX - 0.5));
        }

        // Acceptance criterion
        double acceptanceRatio = distribution(x_new) / distribution(x);
        if (acceptanceRatio >= 1.0 || ((double)rand() / RAND_MAX) < acceptanceRatio) {
            x = x_new;
            if (i >= M) {  // Start counting accepted moves after M steps
                accepted++;
            }
        }
    }

    // Calculate acceptance rate R
    double R = static_cast<double>(accepted) / (N - M);

    // Adjust delta until R is between 0.3 and 0.5
    while (R < 0.3 || R > 0.5) {
        delta *= (R < 0.3) ? k : 1.0 / k;  // Increase or decrease delta
        accepted = 0;
        x = start;

        for (int i = 0; i < N; ++i) {
            std::vector<double> x_new = x;
            for (int j = 0; j < dimensions; ++j) {
                x_new[j] = x[j] + (2.0 * delta * ((double)rand() / RAND_MAX - 0.5));
            }

            double acceptanceRatio = distribution(x_new) / distribution(x);
            if (acceptanceRatio >= 1.0 || ((double)rand() / RAND_MAX) < acceptanceRatio) {
                x = x_new;
                if (i >= M) {
                    accepted++;
                }
            }
        }
        R = static_cast<double>(accepted) / (N - M);
    }

    // Print and return the updated delta
    std::cout << "Adjusted delta: " << delta << std::endl;
    return delta;
}