#include <iostream>
#include <vector>
#include <chrono>

using namespace std;
using namespace chrono;

template <typename T>
T dotProduct(const vector<T>& v1, const vector<T>& v2) {
    T result = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

template <typename T>
vector<T> matVecProduct(const vector<vector<T>>& mat, const vector<T>& vec) {
    vector<T> result(mat.size(), 0);
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[i].size(); ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}

template <typename T>
vector<vector<T>> matMatProduct(const vector<vector<T>>& mat1, const vector<vector<T>>& mat2) {
    vector<vector<T>> result(mat1.size(), vector<T>(mat2[0].size(), 0));
    for (size_t i = 0; i < mat1.size(); ++i) {
        for (size_t j = 0; j < mat2[0].size(); ++j) {
            for (size_t k = 0; k < mat1[0].size(); ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

template <typename T>
void measureTime() {
    size_t n = 1000;
    size_t iterations = 10000;

    // Create random vectors and matrices
    vector<T> vec1(n, 1.0);
    vector<T> vec2(n, 2.0);
    vector<vector<T>> mat1(n, vector<T>(n, 1.0));
    vector<vector<T>> mat2(n, vector<T>(n, 2.0));

    // Dot Product
    auto start = high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        dotProduct(vec1, vec2);
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "C++ Dot Product Duration: " << duration.count() << " microseconds\n";

    // Matrix-Vector Product
    start = high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        matVecProduct(mat1, vec1);
    }
    end = high_resolution_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout << "C++ Matrix-Vector Product Duration: " << duration.count() << " microseconds\n";

    // Matrix-Matrix Product
    start = high_resolution_clock::now();
    for (size_t i = 0; i < iterations; ++i) {
        matMatProduct(mat1, mat2);
    }
    end = high_resolution_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout << "C++ Matrix-Matrix Product Duration: " << duration.count() << " microseconds\n";
}

int main() {
    measureTime<float>(); // For float type
    return 0;
}
