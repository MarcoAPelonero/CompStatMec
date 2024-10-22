#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "metropolis.hpp"
#include "monte_carlo.hpp"
#include "utils.hpp"
#include "distributions.hpp"
#include "progress_bar.hpp"  // Include the custom ProgressBar

// Main function
int main() {
    // Monte Carlo and Metropolis settings
    int dimensions = 2;  // 2D problem (x, y)
    int N_steps = 10000;  // Total steps
    int M_steps = 10;  // Start calculating acceptance rate after M steps
    double initial_delta = 0.1;  // Initial guess for delta
    double k = 1.1;  // Factor for adjusting delta

    // Define the bounds of the map
    int L = 100;
    std::vector<double> lower_bound = {0.0, 0.0};
    std::vector<double> upper_bound = {static_cast<double>(L), static_cast<double>(L)};
    std::vector<double> start = {L / 2.0, L / 2.0};  // Start in the middle of the map

    // Step 1: Evaluate a good delta for Metropolis using utils function
    std::cout << "Evaluating optimal delta for Metropolis..." << std::endl;
    //double optimal_delta = adjustDelta(nile_distribution, start, N_steps, M_steps, initial_delta, k);
    double optimal_delta = 10;
    // Step 2: Initialize Monte Carlo and Metropolis integrators
    MonteCarlo mc(dimensions);
    Metropolis metropolis(optimal_delta, dimensions);

    // Arrays to store integration results at each step
    std::vector<double> mc_results(N_steps);
    std::vector<double> metropolis_results(N_steps);

    // Step 3: Progress bar setup
    ProgressBar mc_bar(N_steps);
    ProgressBar metropolis_bar(N_steps);

    // Step 4: Perform Monte Carlo integration and store results
    std::cout << "Running Monte Carlo integration..." << std::endl;
    for (int i = 0; i < N_steps; ++i) {
        double mc_area = mc.integrate(nile_distribution, lower_bound, upper_bound, i+1);
        mc_results[i] = mc_area;
        mc_bar.update(i + 1);  // Update progress bar
    }
    mc_bar.finish();  // Finish progress bar

    // Step 5: Perform Metropolis integration and store results
    std::cout << "Running Metropolis integration..." << std::endl;
    for (int i = 0; i < N_steps; ++i) {
        double metropolis_area = metropolis.integrate(nile_distribution, start, i+1);
        metropolis_results[i] = metropolis_area;
        metropolis_bar.update(i + 1);  // Update progress bar
    }
    metropolis_bar.finish();  // Finish progress bar

    // Step 6: Calculate the actual area of the Nile
    int nile_width = 10;
    int nile_height = 20;
    double actual_nile_area = nile_width * nile_height;

    // Step 7: Print the final results
    std::cout << "Results after " << N_steps << " steps:" << std::endl;
    std::cout << "Monte Carlo calculated area: " << mc_results[N_steps - 1] << std::endl;
    std::cout << "Metropolis calculated area: " << metropolis_results[N_steps - 1] << std::endl;
    std::cout << "Actual Nile area: " << actual_nile_area << std::endl;

    return 0;
}