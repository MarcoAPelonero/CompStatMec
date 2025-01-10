#include "RDFCalculator.hpp"
#include "boundaryConditions.hpp" // For minimum image convention
#include <cmath>
#include <algorithm>

// Constructor
RDFCalculator::RDFCalculator(int bins, double cutoff)
    : rdfBins(bins), rdfCutoff(cutoff), rdfBinWidth(cutoff / bins), rdfHistogram(bins, 0.0)
{}

// Reset Histogram
void RDFCalculator::resetHistogram() {
    std::fill(rdfHistogram.begin(), rdfHistogram.end(), 0.0);
}

// Direct Computation of RDF
void RDFCalculator::computeDirect(const std::vector<Particle> &particles, double boxLength) {
    resetHistogram();

    for(size_t i = 0; i < particles.size(); ++i) {
        for(size_t j = i + 1; j < particles.size(); ++j) {
            Vec dr = particles[i].getPosition() - particles[j].getPosition();
            // Apply Minimum Image Convention
            double halfBox = boxLength / 2.0;
            if (dr.getX() > halfBox) dr.setX(dr.getX() - boxLength);
            if (dr.getX() < -halfBox) dr.setX(dr.getX() + boxLength);
            if (dr.getY() > halfBox) dr.setY(dr.getY() - boxLength);
            if (dr.getY() < -halfBox) dr.setY(dr.getY() + boxLength);
            if (dr.getZ() > halfBox) dr.setZ(dr.getZ() - boxLength);
            if (dr.getZ() < -halfBox) dr.setZ(dr.getZ() + boxLength);

            double r = dr.length();
            if(r < rdfCutoff) {
                int bin = static_cast<int>(std::floor(r / rdfBinWidth));
                if(bin >=0 && bin < rdfBins) {
                    rdfHistogram[bin] += 2.0; // Each pair counted twice
                }
            }
        }
    }

    // Normalize RDF
    double volume = std::pow(boxLength, 3);
    double density = static_cast<double>(particles.size()) / volume;

    for(int i = 0; i < rdfBins; ++i) {
        double rLower = i * rdfBinWidth;
        double rUpper = rLower + rdfBinWidth;
        double shellVolume = (4.0 / 3.0) * M_PI * (std::pow(rUpper, 3) - std::pow(rLower, 3));
        double idealCount = density * shellVolume * static_cast<double>(particles.size());
        rdfHistogram[i] /= idealCount;
    }
}

// Optimized RDF Computation using Cell Lists (Implementation Omitted for Brevity)
void RDFCalculator::computeOptimized(const std::vector<Particle> &particles, double boxLength) {
    // Implement cell list optimization here
    // This requires additional structures and is omitted for brevity
}

// Compute Thermodynamic Properties from RDF
void RDFCalculator::computeThermodynamics(std::ofstream &file, int numParticles, double boxLength) {
    double temperature = 0.0; // Placeholder: Compute if needed
    double volume = std::pow(boxLength, 3);
    double density = static_cast<double>(numParticles) / volume;

    double dr = rdfBinWidth;
    double energy_per_particle = 0.0;
    double pressure_term = 0.0;
    double k_B = 1.0; // Boltzmann constant

    for(int i = 0; i < rdfBins; ++i) {
        double r = (i + 0.5) * dr;
        double g_r = rdfHistogram[i];

        double r2 = r * r;
        double inv_r2 = 1.0 / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;

        double u_r = 4.0 * (inv_r12 - inv_r6);
        energy_per_particle += g_r * u_r * (r2) * dr;

        double du_dr = 24.0 * (2.0 * inv_r12 - inv_r6) / r;

        pressure_term += g_r * (4.0 * M_PI * std::pow(r, 3) * du_dr) * dr;
    }

    energy_per_particle *= 2.0 * M_PI * density;
    double total_potential_energy = energy_per_particle * numParticles;

    double p = density * k_B * temperature - (2.0 / 3.0) * M_PI * density * density * pressure_term;

    file << temperature << " " << energy_per_particle << " " << total_potential_energy << " " << p << "\n";
}

double RDFCalculator::computeAverageEnergyPerParticle(int numParticles, double boxLength) {
    double volume = std::pow(boxLength, 3);
    double density = static_cast<double>(numParticles) / volume;

    double dr = rdfBinWidth;
    double energy_per_particle = 0.0;
    double pressure_term = 0.0;

    for(int i = 0; i < rdfBins; ++i) {
        double r = (i + 0.5) * dr;
        double g_r = rdfHistogram[i];

        double r2 = r * r;
        double inv_r2 = 1.0 / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;

        double u_r = 4.0 * (inv_r12 - inv_r6);
        energy_per_particle += g_r * u_r * (r2) * dr;

        double du_dr = 24.0 * (2.0 * inv_r12 - inv_r6) / r;

        pressure_term += g_r * (4.0 * M_PI * std::pow(r, 3) * du_dr) * dr;
    }

    energy_per_particle *= 2.0 * M_PI * density;

    return energy_per_particle;
}

void RDFCalculator::print(std::ofstream &file) const {
    for(int i = 0; i < rdfBins; ++i) {
        double r = (i + 0.5) * rdfBinWidth;
        file << r << " " << rdfHistogram[i] << "\n";
    }
    file << "\n";
}

Eigen::RowVectorXf RDFCalculator::getHistogramAsEigenVector() const {
    Eigen::RowVectorXf rdf_vector(rdfBins);
    for (int i = 0; i < rdfBins; ++i) {
        rdf_vector(i) = rdfHistogram[i];
    }
    return rdf_vector;
}