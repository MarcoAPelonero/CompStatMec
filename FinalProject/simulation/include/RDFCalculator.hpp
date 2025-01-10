#ifndef RDF_CALCULATOR_HPP
#define RDF_CALCULATOR_HPP

#include "particle.hpp"
#include "Vec.hpp"

#include <vector>
#include <fstream>
#include <eigen3/Eigen/Dense>

class RDFCalculator {
private:
    int rdfBins;
    double rdfCutoff;
    double rdfBinWidth;
    std::vector<double> rdfHistogram;

public:
    RDFCalculator(int bins, double cutoff);
    void resetHistogram();
    void computeDirect(const std::vector<Particle> &particles, double boxLength);
    void computeOptimized(const std::vector<Particle> &particles, double boxLength);
    void computeThermodynamics(std::ofstream &file, int numParticles, double boxLength);
    double computeAverageEnergyPerParticle(int numParticles, double boxLength);
    void print(std::ofstream &file) const;

    // ML Inference 
    Eigen::RowVectorXf getHistogramAsEigenVector() const;
};

#endif // RDF_CALCULATOR_HPP
