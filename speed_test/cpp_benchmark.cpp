#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <thread>
#include <fstream>

class Vector {
public:
    std::vector<double> data;
    int size;

    Vector(int size) : size(size), data(size) {}

    void randomize() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < size; ++i) {
            data[i] = dis(gen);
        }
    }

    double dot(const Vector& other) const {
        double result = 0.0;
        for (int i = 0; i < size; ++i) {
            result += data[i] * other.data[i];
        }
        return result;
    }

    double dot_optimized(const Vector& other) const {
        double result = 0.0;
        int i = 0;
        int limit = size - (size % 4);

        for (; i < limit; i += 4) {
            result += data[i] * other.data[i];
            result += data[i + 1] * other.data[i + 1];
            result += data[i + 2] * other.data[i + 2];
            result += data[i + 3] * other.data[i + 3];
        }

        for (; i < size; ++i) {
            result += data[i] * other.data[i];
        }

        return result;
    }
};

class Matrix {
public:
    std::vector<std::vector<double>> data;
    int rows, cols;

    Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows, std::vector<double>(cols)) {}

    void randomize() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                data[i][j] = dis(gen);
            }
        }
    }

    Vector multiply(const Vector& vec) const {
        Vector result(rows);
        for (int i = 0; i < rows; ++i) {
            result.data[i] = 0.0;
            for (int j = 0; j < cols; ++j) {
                result.data[i] += data[i][j] * vec.data[j];
            }
        }
        return result;
    }

    Vector multiply_optimized(const Vector& vec) const {
        Vector result(rows);
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            int j = 0;
            int limit = cols - (cols % 4);

            for (; j < limit; j += 4) {
                sum += data[i][j] * vec.data[j];
                sum += data[i][j + 1] * vec.data[j + 1];
                sum += data[i][j + 2] * vec.data[j + 2];
                sum += data[i][j + 3] * vec.data[j + 3];
            }

            for (; j < cols; ++j) {
                sum += data[i][j] * vec.data[j];
            }

            result.data[i] = sum;
        }
        return result;
    }

    Matrix multiply(const Matrix& mat) const {
        Matrix result(rows, mat.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < mat.cols; ++j) {
                result.data[i][j] = 0.0;
                for (int k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * mat.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix multiply_optimized(const Matrix& mat) const {
        Matrix result(rows, mat.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < mat.cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < cols; ++k) {
                    sum += data[i][k] * mat.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        return result;
    }
};

class ProgressBar {
public:
    ProgressBar(int total, int width = 50) : total(total), width(width), progress(0) {}

    void update(int value) {
        progress = value;
        int pos = width * progress / total;
        std::cout << "[";
        for (int i = 0; i < width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0 / total) << " %\r";
        std::cout.flush();
    }

    void finish() {
        update(total);
        std::cout << std::endl;
    }

private:
    int total;
    int width;
    int progress;
};

int main() {
    using namespace std::chrono;
    std::cout << "C++ Benchmark Test" << std::endl;
    const int runs = 100;
    const int n = 1000;
    Vector v1(n), v2(n);
    v1.randomize();
    v2.randomize();
    Matrix M(n, n);
    M.randomize();
    Matrix M2(n, n);
    M2.randomize();

    double scalar_product_time = 0.0;
    double matrix_vector_time = 0.0;
    double matrix_matrix_time = 0.0;

    ProgressBar progressBar(runs);

    for (int i = 0; i < runs; ++i) {
        auto start = high_resolution_clock::now();
        double scalar_product = v1.dot_optimized(v2);
        auto end = high_resolution_clock::now();
        scalar_product_time += duration_cast<microseconds>(end - start).count();

        start = high_resolution_clock::now();
        Vector result_vec = M.multiply_optimized(v1);
        end = high_resolution_clock::now();
        matrix_vector_time += duration_cast<microseconds>(end - start).count();

        start = high_resolution_clock::now();
        Matrix result_mat = M.multiply_optimized(M2);
        end = high_resolution_clock::now();
        matrix_matrix_time += duration_cast<microseconds>(end - start).count();

        progressBar.update(i + 1);
    }

    progressBar.finish();

    std::ofstream file("cpp_benchmark.dat");
    file << "Scalar Product Time: " << scalar_product_time / runs << std::endl;
    file << "Matrix-Vector Multiplication Time: " << matrix_vector_time / runs << std::endl;
    file << "Matrix-Matrix Multiplication Time: " << matrix_matrix_time / runs << std::endl;
    file.close();

    return 0;
}