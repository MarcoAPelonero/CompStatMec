#include "MolecularDynamicsSimulation.hpp"

void MolecularDynamicsSimulation::run() {

    std::ofstream trajectoryFile(outputFileName);
    if (!trajectoryFile.is_open()) {
        std::cerr << "Error: could not open " << outputFileName << " for writing." << std::endl;
        return;
    } 
    
    std::ofstream rdfFile;
    // Print the name of the RDF file 
    std::cout << "RDF file: " << outputRDFFileName << std::endl;
    if (outputRDFFileName != "") {
        
        rdfFile.open(outputRDFFileName);
        if (!rdfFile.is_open()) {
            std::cerr << "Error: could not open " << outputRDFFileName << " for writing." << std::endl;
            return;
        } 
    }

    std::ofstream mlFile;
    if (outputMLFileName != "") {
        mlFile.open(outputMLFileName);
        if (!mlFile.is_open()) {
            std::cerr << "Error: could not open " << outputMLFileName << " for writing." << std::endl;
            return;
        }
    }   

    switch (integrationMethod) {
        case IntegrationMethod::Euler:
            std::cout << "Running Euler" << std::endl;
            runEuler(trajectoryFile);
            break;
        case IntegrationMethod::EulerCromer:
            std::cout << "Running EulerCromer" << std::endl;    
            runEulerCromer(trajectoryFile);
            break;
        case IntegrationMethod::SpeedVerlet:
            std::cout << "Running SpeedVerlet" << std::endl;
            runSpeedVerlet(trajectoryFile);
            break;
        case IntegrationMethod::ThermoEulerCromer:
            std::cout << "Running ThermoEulerCromer" << std::endl;
            runThermoEulerCromer(trajectoryFile);
            break;
        case IntegrationMethod::ThermoSpeedVerlet:
            std::cout << "Running ThermoSpeedVerlet" << std::endl;
            runThermoSpeedVerlet(trajectoryFile);
            break;
        case IntegrationMethod::NPTEulerCromer:
            std::cout << "Running NPTEulerCromer" << std::endl;
            runNPTEulerCromer(trajectoryFile);
            break;
        case IntegrationMethod::NPTSpeedVerlet:
            std::cout << "Running NPTSpeedVerlet" << std::endl;
            runNPTSpeedVerlet(trajectoryFile);
            break;
        case IntegrationMethod::RDFEulerCromer:
            std::cout << "Running RDFEulerCromer" << std::endl;
            runRDFEulerCromer(trajectoryFile, rdfFile);
            break;
        case IntegrationMethod::RDFSpeedVerlet:
            std::cout << "Running RDFSpeedVerlet" << std::endl;
            runRDFSpeedVerlet(trajectoryFile, rdfFile);
            break;
        case IntegrationMethod::MLSpeedVerlet:
            std::cout << "Running MLSpeedVerlet" << std::endl;
            runMLSpeedVerlet(mlFile);
            break;
        default:
            std::cerr << "Error: unknown integration method." << std::endl;
            break;
    }

    if (trajectoryFile.is_open()) trajectoryFile.close();
    if (rdfFile.is_open()) rdfFile.close();
    if (mlFile.is_open()) mlFile.close();
}

void MolecularDynamicsSimulation::runEuler(std::ofstream &trajectoryFile) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectoryFile);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleStepEuler(timeStep);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runEulerCromer(std::ofstream &trajectoryFile) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectoryFile);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleStepEulerCromer(timeStep);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runSpeedVerlet(std::ofstream &trajectoryFile) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectoryFile);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleStepSpeedVerlet(timeStep);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runThermoEulerCromer(std::ofstream &trajectoryFile) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectoryFile);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleStepEulerCromer(timeStep);
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runThermoSpeedVerlet(std::ofstream &trajectoryFile) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectoryFile);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensembleStepSpeedVerlet(timeStep);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    } 
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runNPTEulerCromer(std::ofstream &trajectory) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectory);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleStepEulerCromer(timeStep);
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensemblePressurize(pressure, timeStep, taup);
        ensemble.ensembleSnapshot(trajectory);
        progressBar.update(i);
    }
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runNPTSpeedVerlet(std::ofstream &trajectory) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectory);
    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensemblePressurize(pressure, timeStep, taup);
        ensemble.ensembleStepSpeedVerlet(timeStep);
        ensemble.ensembleSnapshot(trajectory);
        progressBar.update(i);
    }
    progressBar.finish(); 
}

void MolecularDynamicsSimulation::runRDFEulerCromer(std::ofstream &trajectoryFile, std::ofstream &rdfFile) {
    ProgressBar progressBar(totalSteps);
    ensemble.ensembleSnapshot(trajectoryFile);

    auto start_time = std::chrono::high_resolution_clock::now(); // Start timing

    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensembleStepEulerCromer(timeStep);
        ensemble.computeRadialDistributionFunctionDirect();
        rdfFile << "# Step " << i << "\n";
        ensemble.computeThermodynamicsFromRDF(rdfFile);
        ensemble.printRadialDistributionFunction(rdfFile);
        ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish();

    auto end_time = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Simulation completed in " << elapsed.count() << " seconds." << std::endl;
}

void MolecularDynamicsSimulation::runRDFSpeedVerlet(std::ofstream &trajectoryFile, std::ofstream &rdfFile) {
    ProgressBar progressBar(totalSteps);
    // ensemble.ensembleSnapshot(trajectoryFile);

    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensembleStepSpeedVerlet(timeStep);
        ensemble.computeRadialDistributionFunctionDirect();
        rdfFile << "# Step " << i << "\n";
        ensemble.computeThermodynamicsFromRDF(rdfFile);
        ensemble.printRadialDistributionFunction(rdfFile);
        // ensemble.ensembleSnapshot(trajectoryFile);
        progressBar.update(i);
    }
    progressBar.finish();
}

void MolecularDynamicsSimulation::runMLSpeedVerlet(std::ofstream &mlFile) {
    ProgressBar progressBar(totalSteps);

    for (int i = 0; i < totalSteps; ++i) {
        ensemble.ensembleThermalize(temperature, timeStep, taup);
        ensemble.ensembleStepSpeedVerlet(timeStep);
        ensemble.computeRadialDistributionFunctionDirect();

        // 1) Compute thermodynamics from RDF
        auto startRdf = std::chrono::high_resolution_clock::now();
        double rdfEnergy = ensemble.computeEnergyPerParticleFromRDF();
        auto endRdf = std::chrono::high_resolution_clock::now();
        double rdfTime = std::chrono::duration<double>(endRdf - startRdf).count();

        // 2) Compute thermodynamics from ML
        auto startML = std::chrono::high_resolution_clock::now();
        double mlEnergy = ensemble.computeThermodynamicsFromML();
        auto endML = std::chrono::high_resolution_clock::now();
        double mlTime = std::chrono::duration<double>(endML - startML).count();

        // 3) Compute thermodynamics from Python ML
        auto startPyML = std::chrono::high_resolution_clock::now();
        auto [pyEnergy, pyTimeReported] = ensemble.computeThermodynamicsFromPythonML();
        auto endPyML = std::chrono::high_resolution_clock::now();
        double pyTime = std::chrono::duration<double>(endPyML - startPyML).count();

        if (pyEnergy < -1.0 && pyTimeReported < -1.0) {
            std::cerr << "Error in ML prediction" << std::endl;
            return;
        }

        // Write timing and energy values to mlFile
        mlFile << rdfTime << " " << rdfEnergy << " "
               << mlTime << " " << mlEnergy << " "
               << pyTimeReported << " " << pyEnergy << std::endl;

        progressBar.update(i);
    }
    progressBar.finish();
}