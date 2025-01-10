// include/rdf_density_model.hpp

#pragma once

#include <eigen3/Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

class RDFDensityModel {
public:
    // Model Dimensions
    static const int RDF_SIZE = 1000;
    static const int HIDDEN_DIM = 128;
    static const int DENSITY_DIM = 1;
    static const int DENSITY_HIDDEN_DIM = HIDDEN_DIM / 2; // 64
    static const int MERGED_DIM = HIDDEN_DIM + DENSITY_HIDDEN_DIM; // 192
    static const int OUTPUT_DIM = 1;

    // Constructor
    RDFDensityModel();

    // Load weights
    void load_weights(const std::string& folder_path);

    // Predict energy
    float predict_energy(const Eigen::RowVectorXf& rdf_input, const Eigen::RowVectorXf& density_input);

private:
    // Helper methods
    Eigen::MatrixXf loadCSV(const std::string& csvPath, int rows, int cols);
    Eigen::RowVectorXf relu(const Eigen::RowVectorXf& inputRow);

    // Weight and bias matrices
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
};
