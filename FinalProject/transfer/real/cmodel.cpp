// rdf_density_model.cpp

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * @brief Load a CSV file into an Eigen::MatrixXf.
 *
 * @param csvPath Path to the CSV file.
 * @param rows Number of rows in the CSV.
 * @param cols Number of columns in the CSV.
 * @return Eigen::MatrixXf Loaded matrix.
 */
Eigen::MatrixXf loadCSV(const std::string& csvPath, int rows, int cols) {
    std::ifstream file(csvPath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << csvPath << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<float> values;
    values.reserve(rows * cols);

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stof(cell));
        }
    }
    file.close();

    if (values.size() != static_cast<size_t>(rows * cols)) {
        std::cerr << "Error: Expected " << rows * cols << " values, but got " << values.size() << " in " << csvPath << std::endl;
        exit(1);
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

/**
 * @brief Apply ReLU activation to an Eigen::RowVectorXf.
 *
 * @param inputRow Input vector.
 * @return Eigen::RowVectorXf Activated vector.
 */
Eigen::RowVectorXf relu(const Eigen::RowVectorXf& inputRow) {
    return inputRow.cwiseMax(0.0f);
}

// Define model dimensions
const int RDF_SIZE = 1000;
const int HIDDEN_DIM = 128;
const int DENSITY_DIM = 1;
const int DENSITY_HIDDEN_DIM = HIDDEN_DIM / 2; // 64
const int MERGED_DIM = HIDDEN_DIM + DENSITY_HIDDEN_DIM; // 192
const int OUTPUT_DIM = 1;

// Weight and bias matrices (using float)
Eigen::MatrixXf rdf_branch_0_weight;    // (128, 1000)
Eigen::VectorXf rdf_branch_0_bias;      // (128)

Eigen::MatrixXf rdf_branch_2_weight;    // (128, 128)
Eigen::VectorXf rdf_branch_2_bias;      // (128)

Eigen::MatrixXf density_branch_0_weight; // (64, 1)
Eigen::VectorXf density_branch_0_bias;   // (64)

Eigen::MatrixXf density_branch_2_weight; // (64, 64)
Eigen::VectorXf density_branch_2_bias;   // (64)

Eigen::MatrixXf merged_head_0_weight;   // (128, 192)
Eigen::VectorXf merged_head_0_bias;     // (128)

Eigen::MatrixXf merged_head_2_weight;   // (1, 128)
Eigen::VectorXf merged_head_2_bias;     // (1)

/**
 * @brief Load all weights and biases from CSV files within the specified folder.
 *
 * @param folder_path Path to the 'layers' folder containing CSV files.
 */
void load_weights(const std::string& folder_path) {
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

/**
 * @brief Forward pass for the RDF branch.
 *
 * @param x_rdf Input RDF vector (1 x 1000).
 * @return Eigen::RowVectorXf Output after RDF branch (1 x 128).
 */
Eigen::RowVectorXf forward_rdf(const Eigen::RowVectorXf& x_rdf) {
    // Layer 0: Linear + ReLU
    Eigen::RowVectorXf out = (x_rdf * rdf_branch_0_weight.transpose()) + rdf_branch_0_bias.transpose();
    out = relu(out); // (1 x 128)

    // Layer 2: Linear + ReLU
    out = (out * rdf_branch_2_weight.transpose()) + rdf_branch_2_bias.transpose();
    out = relu(out); // (1 x 128)

    return out; // (1 x 128)
}

/**
 * @brief Forward pass for the Density branch.
 *
 * @param x_density Input density (1 x 1).
 * @return Eigen::RowVectorXf Output after Density branch (1 x 64).
 */
Eigen::RowVectorXf forward_density(const Eigen::RowVectorXf& x_density) {
    // Layer 0: Linear + ReLU
    Eigen::RowVectorXf out = (x_density * density_branch_0_weight.transpose()) + density_branch_0_bias.transpose();
    out = relu(out); // (1 x 64)

    // Layer 2: Linear + ReLU
    out = (out * density_branch_2_weight.transpose()) + density_branch_2_bias.transpose();
    out = relu(out); // (1 x 64)

    return out; // (1 x 64)
}

/**
 * @brief Forward pass for the Merged head.
 *
 * @param merged_input Concatenated vector from RDF and Density branches (1 x 192).
 * @return Eigen::RowVectorXf Output energy prediction (1 x 1).
 */
Eigen::RowVectorXf forward_merged(const Eigen::RowVectorXf& merged_input) {
    // Layer 0: Linear + ReLU
    Eigen::RowVectorXf out = (merged_input * merged_head_0_weight.transpose()) + merged_head_0_bias.transpose();
    out = relu(out); // (1 x 128)

    // Layer 2: Linear (no activation)
    out = (out * merged_head_2_weight.transpose()) + merged_head_2_bias.transpose();
    // No activation on the final layer

    return out; // (1 x 1)
}

int main() {
    // 1. Define the path to the 'layers' folder
    std::string layers_dir = "layers";

    // 2. Load all weights and biases from the 'layers' folder
    load_weights(layers_dir);

    // 3. Define input data
    // Example inputs:
    // - x_rdf: [0.5, 0.5, ..., 0.5] (1000 elements)
    // - x_density: [1.0]

    Eigen::RowVectorXf x_rdf(RDF_SIZE);
    x_rdf.setConstant(0.5f); // All entries set to 0.5

    Eigen::RowVectorXf x_density(DENSITY_DIM);
    x_density << 1.0f; // Single density value set to 1.0

    // 4. Forward pass through RDF branch
    Eigen::RowVectorXf out_rdf = forward_rdf(x_rdf);

    // 5. Forward pass through Density branch
    Eigen::RowVectorXf out_density = forward_density(x_density);

    // 6. Concatenate outputs
    Eigen::RowVectorXf merged(MERGED_DIM);
    merged << out_rdf, out_density;

    // 7. Forward pass through Merged head
    Eigen::RowVectorXf energy = forward_merged(merged);

    // 8. Print the output
    std::cout << "C++ forward pass output (Energy prediction): " << energy(0) << std::endl;

    return 0;
}
