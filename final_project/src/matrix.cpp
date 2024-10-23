#include "matrix.hpp"

Matrix::Matrix() { for (int i = 0; i < N1; ++i) for (int j = 0; j < N2; ++j) {M[i][j] = 0;} }

Matrix::Matrix(std::initializer_list<ntype> rows, std::initializer_list<ntype> cols) {
    if (rows.size() != N1 || cols.size() != N2)
        throw std::invalid_argument("Initializer list size does not match matrix dimensions.");
    int r = 0;
    for (auto row : rows) {
        int c = 0;
        for (auto col : cols) {
            M[r][c] = row * col;
            c++;
        }
        r++;
    }
}

Matrix::Matrix(std::initializer_list<ntype> diag) {
    if (diag.size() != N1 || N1 != N2)
        throw std::invalid_argument("Diagonal initializer requires a square matrix.");
    for (int i = 0; i < N1; i++)
        for (int j = 0; j < N2; j++) 
            M[i][j] = 0;
    int c = 0, d = 0;
    for (auto el : diag) 
        M[d++][c++] = el;
}

Matrix::~Matrix() {}

void Matrix::show(std::string name) const {
    std::cout << name << std::endl;
    for (int i = 0; i < N1; ++i) {
        std::cout << "| ";
        for (int j = 0; j < N2; ++j) {
            std::cout << M[i][j];
            std::cout << " ";
        }
        std::cout << "|" << std::endl;
    }
}

ntype Matrix::get(int i, int j) const {
    return M[i][j];
}

ntype Matrix::set(int i, int j, ntype val) {
    return M[i][j] = val;
}

Matrix Matrix::diagonalize(int &swaps) const {
    if (N1 != N2) {
        throw std::invalid_argument("Matrix must be square.");
    }

    // Copy the original matrix so we may not modify the original element but cretae a new element
    Matrix U(*this);  
    swaps = 0;
    for (int i = 0; i < N1; i++) {
        // Find the pivot element
        if (U(i, i) == 0) {
            // Find a row below to swap with
            bool swapped = false;
            for (int j = i + 1; j < N1; j++) {
                if (U(j, i) != 0) {
                    std::swap(U.M[i], U.M[j]);
                    swaps++;
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                return U;  // Singular matrix, return early
            }
        }

        // Gaussian elimination to create upper triangular matrix
        for (int j = i + 1; j < N1; j++) {
            double factor = U(j, i) / U(i, i);
            for (int k = i; k < N1; k++) {
                U(j, k) -= factor * U(i, k);
            }
        }
    }

    return U;
}

double Matrix::determinant() const {
    int swaps = 0;
    Matrix U = diagonalize(swaps);

    // Calculate the determinant by multiplying the diagonal elements
    double det = 1;
    for (int i = 0; i < N1; i++) {
        det *= U(i, i);
    }

    // Adjust for row swaps
    if (swaps % 2 != 0) {
        det = -det;
    }

    return det;
}



Matrix& Matrix::operator=(const Matrix& M2) {
    for (int i = 0; i < N1; ++i) 
        for (int j = 0; j < N2; ++j)
                M[i][j] = M2.M[i][j];
    
    return *this;
}

Matrix Matrix::operator+(const Matrix& M2) {
    Matrix M3;
    for (int i = 0; i < N1; ++i) 
        for (int j = 0; j < N2; ++j)
            M3.M[i][j] = M[i][j] + M2.M[i][j];
    
    return M3;
}

Matrix Matrix::operator-(const Matrix& M2) {
    Matrix M3;
    for (int i = 0; i < N1; ++i) 
        for (int j = 0; j < N2; ++j)
            M3.M[i][j] = M[i][j] - M2.M[i][j];
    
    return M3;
}

Matrix Matrix::operator*(const Matrix& M2) {
    Matrix M3;
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            ntype el = 0;
            for (int k = 0; k < N1; ++k) {
                el += M[i][k] * M2.M[k][j];
            }
            M3.M[i][j] = el;
        }
    }
    return M3;
}


ntype& Matrix::operator()(int i, int j) {
        return M[i][j];
    }