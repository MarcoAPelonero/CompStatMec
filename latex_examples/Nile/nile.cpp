#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>

const int gridSize = 100;
const double L = 1.0;
const double depthOutside = 0.0;
const double depthMax = 10.0;

double depthFunction(int x, int y) {
    return depthOutside + (depthMax * std::sin(0.1 * x) * std::cos(0.1 * y));
}

std::vector<std::vector<double>> depth(gridSize, std::vector<double>(gridSize, depthOutside));

void generateNilePath() {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, gridSize - 1);
    int x = gridSize / 2;
    for (int y = 0; y < gridSize; y++) {
        x += distribution(generator) % 3 - 1;
        x = std::max(0, std::min(x, gridSize - 1));
        depth[y][x] = depthFunction(x, y);
    }
}

void monteCarloEstimate(int samples, std::vector<double>& mc_depths) {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, gridSize - 1);
    double sumDepth = 0.0;

    for (int i = 1; i <= samples; i++) {
        int x = distribution(generator);
        int y = distribution(generator);
        sumDepth += depth[y][x];
        mc_depths.push_back(sumDepth / i);  
    }
}

void metropolisEstimate(int steps, std::vector<double>& mcm_depths) {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, gridSize - 1);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    int x = gridSize / 2, y = gridSize / 2;
    double sumDepth = depth[y][x];

    for (int i = 1; i <= steps; i++) {
        int newX = x + distribution(generator) % 3 - 1;
        int newY = y + distribution(generator) % 3 - 1;
        newX = std::max(0, std::min(newX, gridSize - 1));
        newY = std::max(0, std::min(newY, gridSize - 1));

        if (uniform(generator) < std::min(1.0, depth[newY][newX] / depth[y][x])) {
            x = newX;
            y = newY;
            sumDepth += depth[y][x];
        }
        mcm_depths.push_back(sumDepth / i);  
    }
}

int main() {
    generateNilePath();

    int samples = 100000;  // Number of samples for MC
    int steps = 100000;    // Number of steps for MCM

    std::vector<double> mc_depths;
    std::vector<double> mcm_depths;

    monteCarloEstimate(samples, mc_depths);
    metropolisEstimate(steps, mcm_depths);

    // Save results to file
    std::ofstream outfile("nile_depths.txt");
    for (int i = 0; i < samples; i++) {
        outfile << mc_depths[i] << "," << mcm_depths[i] << "\n";
    }
    outfile.close();

    std::cout << "Data saved to nile_depths.txt" << std::endl;
    return 0;
}