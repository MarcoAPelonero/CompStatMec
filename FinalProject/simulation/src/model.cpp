// src/rdf_density_model.cpp

#include "model.hpp"

RDFDensityModel::RDFDensityModel() {
    // Initialization can be done here if necessary
}

Eigen::MatrixXf RDFDensityModel::loadCSV(const std::string& csvPath, int rows, int cols) {
    std::ifstream file(csvPath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << csvPath << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<float> values;
    values.reserve(rows * cols);

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            try {
                values.push_back(std::stof(cell));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid number '" << cell << "' in file " << csvPath << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    file.close();

    if (values.size() != static_cast<size_t>(rows * cols)) {
        std::cerr << "Error: Expected " << rows * cols << " values, but got " << values.size()
                  << " in " << csvPath << std::endl;
        exit(EXIT_FAILURE);
    }

    Eigen::MatrixXf mat(rows, cols);
    int idx = 0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            mat(r, c) = values[idx++];
        }
    }
    return mat;
}

Eigen::RowVectorXf RDFDensityModel::relu(const Eigen::RowVectorXf& inputRow) {
    return inputRow.cwiseMax(0.0f);
}

void RDFDensityModel::load_weights(const std::string& folder_path) {
    std::cout << "Loading weights and biases from '" << folder_path << "' directory...\n";

    // RDF Branch
    rdf_branch_0_weight = loadCSV(folder_path + "/rdf_branch.0_weight.csv", 128, 1000);
    Eigen::MatrixXf temp_bias = loadCSV(folder_path + "/rdf_branch.0_bias.csv", 1, 128);
    rdf_branch_0_bias = temp_bias.row(0).transpose(); // (128)

    rdf_branch_2_weight = loadCSV(folder_path + "/rdf_branch.2_weight.csv", 128, 128);
    temp_bias = loadCSV(folder_path + "/rdf_branch.2_bias.csv", 1, 128);
    rdf_branch_2_bias = temp_bias.row(0).transpose(); // (128)

    // Density Branch
    density_branch_0_weight = loadCSV(folder_path + "/density_branch.0_weight.csv", 64, 1);
    temp_bias = loadCSV(folder_path + "/density_branch.0_bias.csv", 1, 64);
    density_branch_0_bias = temp_bias.row(0).transpose(); // (64)

    density_branch_2_weight = loadCSV(folder_path + "/density_branch.2_weight.csv", 64, 64);
    temp_bias = loadCSV(folder_path + "/density_branch.2_bias.csv", 1, 64);
    density_branch_2_bias = temp_bias.row(0).transpose(); // (64)

    // Merged Head
    merged_head_0_weight = loadCSV(folder_path + "/merged_head.0_weight.csv", 128, 192);
    temp_bias = loadCSV(folder_path + "/merged_head.0_bias.csv", 1, 128);
    merged_head_0_bias = temp_bias.row(0).transpose(); // (128)

    merged_head_2_weight = loadCSV(folder_path + "/merged_head.2_weight.csv", 1, 128);
    temp_bias = loadCSV(folder_path + "/merged_head.2_bias.csv", 1, 1);
    merged_head_2_bias = temp_bias.row(0).transpose(); // (1)

    std::cout << "All weights and biases loaded successfully.\n";
}

float RDFDensityModel::predict_energy(const Eigen::RowVectorXf& rdf_input, const Eigen::RowVectorXf& density_input) {
    // Forward pass through RDF branch
    Eigen::RowVectorXf out_rdf = (rdf_input * rdf_branch_0_weight.transpose()) + rdf_branch_0_bias.transpose();
    out_rdf = relu(out_rdf); // (1 x 128)

    out_rdf = (out_rdf * rdf_branch_2_weight.transpose()) + rdf_branch_2_bias.transpose();
    out_rdf = relu(out_rdf); // (1 x 128)

    // Forward pass through Density branch
    Eigen::RowVectorXf out_density = (density_input * density_branch_0_weight.transpose()) + density_branch_0_bias.transpose();
    out_density = relu(out_density); // (1 x 64)

    out_density = (out_density * density_branch_2_weight.transpose()) + density_branch_2_bias.transpose();
    out_density = relu(out_density); // (1 x 64)

    // Concatenate RDF and Density outputs
    Eigen::RowVectorXf merged(RDFDensityModel::MERGED_DIM);
    merged << out_rdf, out_density; // (1 x 192)

    // Forward pass through Merged head
    Eigen::RowVectorXf out_merged = (merged * merged_head_0_weight.transpose()) + merged_head_0_bias.transpose();
    out_merged = relu(out_merged); // (1 x 128)

    out_merged = (out_merged * merged_head_2_weight.transpose()) + merged_head_2_bias.transpose();
    // No activation on the final layer

    // Assuming output is (1 x 1)
    return out_merged(0);
}
