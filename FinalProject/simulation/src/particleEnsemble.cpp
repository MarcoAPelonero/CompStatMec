// src/ParticleEnsemble.cpp

#include "ParticleEnsemble.hpp"
#include <iostream>
#include <cmath>
// Constructor
ParticleEnsemble::ParticleEnsemble(int N, double L)
    : numParticles(N),
      boxLength(L),
      rng(42), // Fixed seed for reproducibility
      gaussian(0.0, 1.0),
      totalEnergy(0.0),
      totalKineticEnergy(0.0),
      totalVirialEnergy(0.0),
      totalPotentialEnergy(0.0),
      forceCalculator(),
      thermodynamics(),
      rdfCalculator(1000, L/2),
      boundaryConditions(L),
      mlModel() // Initialize ML model
{
    initializeParticles();
    // Load ML model weights
    mlModel.load_weights("layers"); // Assuming 'layers' is in the current working directory
}

// Initialize Particles in a Cubic Lattice with Gaussian Velocities
void ParticleEnsemble::initializeParticles() {
    int numPerSide = std::ceil(std::cbrt(numParticles));
    double spacing = boxLength / numPerSide;

    int particleCount = 0;
    for (int x = 0; x < numPerSide && particleCount < numParticles; ++x) {
        for (int y = 0; y < numPerSide && particleCount < numParticles; ++y) {
            for (int z = 0; z < numPerSide && particleCount < numParticles; ++z) {
                Vec position(x * spacing, y * spacing, z * spacing);
                Vec velocity(gaussian(rng), gaussian(rng), gaussian(rng));
                Particle particle(1.0, position, velocity); 
                particles.push_back(particle);
                ++particleCount;
            }
        }
    }

    // Compute Initial Energies and Forces
    forceCalculator.computePotentials(particles, boxLength);
    thermodynamics.computeKinetics(particles);
    thermodynamics.computeVirial(particles, boxLength, forceCalculator);
    thermodynamics.computeEnergies(totalEnergy, totalKineticEnergy, totalVirialEnergy, totalPotentialEnergy, particles);
}

// Accessors and Mutators
std::vector<Particle>& ParticleEnsemble::getParticles() { return particles; }
const std::vector<Particle>& ParticleEnsemble::getParticles() const { return particles; }
double ParticleEnsemble::getBoxLength() const { return boxLength; }
double ParticleEnsemble::getTotalEnergy() const { return totalEnergy; }
double ParticleEnsemble::getTotalKineticEnergy() const { return totalKineticEnergy; }
double ParticleEnsemble::getTotalPotentialEnergy() const { return totalPotentialEnergy; }

void ParticleEnsemble::setParticles(const std::vector<Particle> &p) { particles = p; }
void ParticleEnsemble::setBoxLength(double l) { boxLength = l; }
void ParticleEnsemble::setTotalEnergy(double e) { totalEnergy = e; }
void ParticleEnsemble::setTotalKineticEnergy(double e) { totalKineticEnergy = e; }
void ParticleEnsemble::setTotalPotentialEnergy(double e) { totalPotentialEnergy = e; }
void ParticleEnsemble::addParticle(const Particle &p) { particles.push_back(p); }

// Display Particle Information
void ParticleEnsemble::showParticles() const {
    for (size_t i = 0; i < particles.size(); ++i) {
        std::cout << "Particle " << i + 1 << ":\n";
        particles[i].show("\t");
    }
}

// Simulation Control Methods
void ParticleEnsemble::ensembleThermalize(double targetTemperature, double dt, double tauT) {
    thermodynamics.thermalize(particles, targetTemperature, dt, tauT);
}

void ParticleEnsemble::ensemblePressurize(double targetPressure, double dt, double tauP) {
    thermodynamics.pressurize(particles, targetPressure, dt, tauP, boxLength, boundaryConditions);
}

// Integration Steps
void ParticleEnsemble::ensembleStepEuler(double dt) {
    EulerIntegrator integrator;
    integrator.step(particles, forceCalculator, boundaryConditions, dt, boxLength);

    // Update Energies
    forceCalculator.computePotentials(particles, boxLength);
    thermodynamics.computeKinetics(particles);
    thermodynamics.computeVirial(particles, boxLength, forceCalculator);
    thermodynamics.computeEnergies(totalEnergy, totalKineticEnergy, totalVirialEnergy, totalPotentialEnergy, particles);
}

void ParticleEnsemble::ensembleStepEulerCromer(double dt) {
    EulerCromerIntegrator integrator;
    integrator.step(particles, forceCalculator, boundaryConditions, dt, boxLength);

    // Update Energies
    forceCalculator.computePotentials(particles, boxLength);
    thermodynamics.computeKinetics(particles);
    thermodynamics.computeVirial(particles, boxLength, forceCalculator);
    thermodynamics.computeEnergies(totalEnergy, totalKineticEnergy, totalVirialEnergy, totalPotentialEnergy, particles);
}

void ParticleEnsemble::ensembleStepSpeedVerlet(double dt) {
    SpeedVerletIntegrator integrator;
    integrator.step(particles, forceCalculator, boundaryConditions, dt, boxLength);

    // Update Energies
    forceCalculator.computePotentials(particles, boxLength);
    thermodynamics.computeKinetics(particles);
    thermodynamics.computeVirial(particles, boxLength, forceCalculator);
    thermodynamics.computeEnergies(totalEnergy, totalKineticEnergy, totalVirialEnergy, totalPotentialEnergy, particles);
}

// Data Output
void ParticleEnsemble::ensembleSnapshot(std::ofstream &file) {
    file << boxLength << "\n";
    for (const auto &particle : particles) {
        file << particle.getPosition() << " "
             << particle.getVelocity() << " "
             << particle.getPotentialEnergy() << " "
             << particle.getKineticEnergy() << " "
             << particle.getTotalEnergy() << "\n";
    }
    file << "\n"; // Separator for snapshots
}

// RDF Computations
void ParticleEnsemble::resetRDFHistogram() {
    rdfCalculator.resetHistogram();
}

void ParticleEnsemble::computeRadialDistributionFunctionDirect() {
    rdfCalculator.computeDirect(particles, boxLength);
}

void ParticleEnsemble::computeRadialDistributionFunctionOptimized() {
    rdfCalculator.computeOptimized(particles, boxLength);
}

void ParticleEnsemble::computeThermodynamicsFromRDF(std::ofstream &file) {
    rdfCalculator.computeThermodynamics(file, numParticles, boxLength);
}

double ParticleEnsemble::computeEnergyPerParticleFromRDF() {
    return rdfCalculator.computeAverageEnergyPerParticle(numParticles, boxLength);
}

double ParticleEnsemble::computeThermodynamicsFromML() {
    // Extract RDF data
    Eigen::RowVectorXf rdf_data = rdfCalculator.getHistogramAsEigenVector(); // Corrected call

    // Prepare density input
    // Assuming density = N / V, where V = boxLength^3 for a cubic box
    double volume = std::pow(boxLength, 3);
    double density = static_cast<double>(numParticles) / volume;
    Eigen::RowVectorXf density_input(RDFDensityModel::DENSITY_DIM);
    density_input << static_cast<float>(density); // Cast to float if necessary

    // Perform ML prediction
    float predicted_energy = mlModel.predict_energy(rdf_data, density_input);

    // Update thermodynamic properties based on ML prediction;
    // Optionally, update other energies if the ML model provides them
    // For example:
    // totalKineticEnergy = ...;
    // totalPotentialEnergy = ...;
    // totalVirialEnergy = ...;

    // Log or output the predicted energy as needed
    // std::cout << "ML Predicted Total Energy: " << totalEnergy << std::endl;
    return predicted_energy;
}



void ParticleEnsemble::printRadialDistributionFunction(std::ofstream &file) {
    rdfCalculator.print(file);
}

using json = nlohmann::json;

std::pair<double, double> ParticleEnsemble::computeThermodynamicsFromPythonML() {
    // 1. Extract RDF data (ensure it has exactly 1000 elements)
    Eigen::RowVectorXf rdf_data = rdfCalculator.getHistogramAsEigenVector();
    
    // 2. Calculate density
    double volume = std::pow(boxLength, 3);
    double density = static_cast<double>(numParticles) / volume;
    
    // 3. Create JSON input
    json input_json;
    input_json["rdf"] = std::vector<float>(rdf_data.data(), rdf_data.data() + rdf_data.size());
    input_json["density"] = static_cast<float>(density);
    std::string input_str = input_json.dump();
    
    // 4. Define paths (use absolute paths if necessary)
    std::string python_script = "scripts/predict_energy.py"; // Update if using absolute path
    std::string model_path = "layers/best_model.pth";        // Update if using absolute path
    
    // 5. Temporary file names
    std::string temp_input_file = "temp_input.json";
    std::string temp_output_file = "temp_output.txt";
    
    // 6. Write JSON to temp_input.json
    std::ofstream input_file(temp_input_file);
    if (!input_file) {
        std::cerr << "Error: Could not open " << temp_input_file << " for writing." << std::endl;
        return std::make_pair(-1.0,-1.0);
    }
    input_file << input_str;
    input_file.close();
    
    // 7. Construct and execute the command
    std::string command = "python " + python_script + " " + model_path + " < " + temp_input_file + " > " + temp_output_file;
    int ret = system(command.c_str());
    if (ret != 0) {
        std::cerr << "Error: Python script failed with return code " << ret << "." << std::endl;
        remove(temp_input_file.c_str());
        remove(temp_output_file.c_str());
        return std::make_pair(-1.0,-1.0);
    }
    
    // 8. Read predicted energy from temp_output.txt
    std::ifstream output_file(temp_output_file);
    if (!output_file) {
        std::cerr << "Error: Could not open " << temp_output_file << " for reading." << std::endl;
        remove(temp_input_file.c_str());
        remove(temp_output_file.c_str());
        return std::make_pair(-1.0,-1.0);
    }
    
    float predicted_energy;
    double inference_time;
    output_file >> predicted_energy >> inference_time;
    if (output_file.fail()) {
        std::cerr << "Error: Failed to read prediction from " << temp_output_file << "." << std::endl;
        remove(temp_input_file.c_str());
        remove(temp_output_file.c_str());
        return std::make_pair(-1.0,-1.0);
    }
    output_file.close();
    
    // 9. Update thermodynamic properties
    totalEnergy = static_cast<double>(predicted_energy);
    
    // 10. Clean up temporary files
    remove(temp_input_file.c_str());
    remove(temp_output_file.c_str());
    
    // Return energy and timing
    return std::make_pair(totalEnergy, inference_time);
}