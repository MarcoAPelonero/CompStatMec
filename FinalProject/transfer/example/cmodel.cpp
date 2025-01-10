#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Helper function: load a CSV file (rows x cols) into an Eigen::MatrixXd.
 *
 * We assume the CSV has exactly 'rows' lines and each line has 'cols' values.
 */
Eigen::MatrixXd loadCSV(const std::string& csvPath, int rows, int cols) {
    std::ifstream file(csvPath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << csvPath << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<double> values;
    values.reserve(rows * cols); // we know how many total values to expect

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
    }
    file.close();

    // Now place the values into an Eigen::MatrixXd
    Eigen::MatrixXd mat(rows, cols);
    int idx = 0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            mat(r, c) = values[idx++];
        }
    }
    return mat;
}

// MLP dimensions must match the Python code
static const int in_dim     = 4;
static const int hidden_dim = 8;
static const int out_dim    = 2;

// We'll store the weights and biases in global variables for convenience
Eigen::MatrixXd fc1_weight;  // shape: (8,4)
Eigen::VectorXd fc1_bias;    // shape: (8)
Eigen::MatrixXd fc2_weight;  // shape: (2,8)
Eigen::VectorXd fc2_bias;    // shape: (2)

/**
 * @brief A custom "model summary" function for this specific 2-layer MLP.
 *
 * We'll manually print out:
 *   - layer name
 *   - input shape
 *   - output shape
 *   - number of parameters
 *
 * It's not dynamic like torchsummary, but you can expand it for bigger models.
 */
void printModelSummary() {
    std::cout << "===== Custom Model Summary (C++) =====\n\n";

    // Layer 1 (FC1)
    {
        // fc1_weight shape: (hidden_dim, in_dim) => (8,4)
        // fc1_bias shape:   (hidden_dim)         => (8)
        const int inShape  = in_dim;
        const int outShape = hidden_dim;
        const int nParams  = (in_dim * hidden_dim) + hidden_dim;
        std::cout << "Layer (FC1)      | " 
                  << "Input: [" << inShape << "]"
                  << "  => Output: [" << outShape << "]"
                  << "  #Params: " << nParams << "\n";
    }

    // ReLU (no parameters, just an activation)
    {
        std::cout << "Layer (ReLU)     | "
                  << "No parameters (activation)\n";
    }

    // Layer 2 (FC2)
    {
        // fc2_weight shape: (out_dim, hidden_dim) => (2,8)
        // fc2_bias shape:   (out_dim)             => (2)
        const int inShape  = hidden_dim;
        const int outShape = out_dim;
        const int nParams  = (hidden_dim * out_dim) + out_dim;
        std::cout << "Layer (FC2)      | "
                  << "Input: [" << inShape << "]"
                  << "  => Output: [" << outShape << "]"
                  << "  #Params: " << nParams << "\n";
    }

    // Totals
    int totalParams = (in_dim * hidden_dim + hidden_dim) + (hidden_dim * out_dim + out_dim);
    std::cout << "\nTotal params: " << totalParams << std::endl << std::endl;
}

/**
 * @brief A simple ReLU function for an Eigen row vector: ReLU(x) = max(0, x)
 */
Eigen::RowVectorXd relu(const Eigen::RowVectorXd& inputRow) {
    Eigen::RowVectorXd result = inputRow;
    for (int i = 0; i < result.size(); ++i) {
        if (result(i) < 0.0) {
            result(i) = 0.0;
        }
    }
    return result;
}

/**
 * @brief Forward pass for a single sample x of shape (1, in_dim).
 *        Returns a row vector of shape (1, out_dim).
 */
Eigen::RowVectorXd forward(const Eigen::RowVectorXd& x) {
    // FC1: x * W1^T + b1
    // x: (1,4), W1: (8,4), so W1^T: (4,8).  But we do x * W1^T => (1,8).
    Eigen::RowVectorXd out1 = x * fc1_weight.transpose(); // shape (1,8)
    // Add bias (8) => need to treat bias as row vector (1,8)
    out1 += fc1_bias.transpose();

    // ReLU
    out1 = relu(out1); // shape (1,8)

    // FC2: out1 * W2^T + b2
    // out1: (1,8), W2: (2,8), W2^T: (8,2). So out2 => (1,2).
    Eigen::RowVectorXd out2 = out1 * fc2_weight.transpose(); // shape (1,2)
    // Add bias (2) => treat as row vector (1,2)
    out2 += fc2_bias.transpose();

    return out2; // shape (1,2)
}

int main() {
    // 1) Load CSVs
    // fc1_weight: shape (8,4)
    fc1_weight = loadCSV("fc1_weight.csv", hidden_dim, in_dim);

    // fc1_bias: shape (1,8), but we want a vector of size (8)
    {
        Eigen::MatrixXd temp = loadCSV("fc1_bias.csv", 1, hidden_dim);
        fc1_bias = temp.row(0).transpose(); // shape(8)
    }

    // fc2_weight: shape (2,8)
    fc2_weight = loadCSV("fc2_weight.csv", out_dim, hidden_dim);

    // fc2_bias: shape (1,2), we want a vector of size (2)
    {
        Eigen::MatrixXd temp = loadCSV("fc2_bias.csv", 1, out_dim);
        fc2_bias = temp.row(0).transpose(); // shape(2)
    }

    // 2) Print a "model summary" in C++
    printModelSummary();

    // 3) Let’s do the same input as in Python: x = [0.5, -1.0, 2.0, 0.0]
    Eigen::RowVectorXd x(in_dim);
    x << 0.5, -1.0, 2.0, 0.0;

    // 4) Run forward pass
    Eigen::RowVectorXd y = forward(x);

    // 5) Print output
    std::cout << "C++ forward pass output for x=[0.5, -1.0, 2.0, 0.0]: " 
              << y << std::endl;

    return 0;
}
