#include "tensors.hpp"
#include "config.hpp"
#include <iostream>
#include <initializer_list>
#include <string>

Scalar::Scalar() : x(0) {}

Scalar::Scalar(ntype num) : x(num) {}

Scalar::~Scalar() {} 

ntype Scalar::get() const {
    return x;
}

Vector::Vector() {
    for (int i = 0; i < N1; ++i) {
        v[i] = 0;
    }
}

Vector::Vector(ntype val) {
    for (int i = 0; i < N1; ++i) {
        v[i] = val;
    }
}

Vector::Vector(std::initializer_list<ntype> list) {
    int c = 0;
    for (auto el : list) {
        if (c < N1) {
            v[c++] = el;
        } else {
            break;
        }
    }
    for (; c < N1; ++c) {
        v[c] = 0.0;
    }
}

Vector::~Vector() {} 

void Vector::show(std::string name) const {
    std::cout << name << "(";
    for (int i = 0; i < N1; ++i) {
        std::cout << v[i];
        if (i < N1 - 1) 
            std::cout << ",";
    }
    std::cout << ")\n";
}

ntype Vector::get(int i) const {
    return v[i];
}

ntype Vector::set(int i, ntype val) {
    return v[i] = val;
}

double Vector::modulus() {
    return (*this) * (*this);
}

Vector& Vector::operator=(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        v[i] = v2.v[i];
    }
    return *this;
}

Vector Vector::operator+(const Vector& v2) {
    Vector v3;
    for (int i = 0; i < N1; ++i) {
        v3.v[i] = v[i] + v2.v[i];
    }
    return v3;
}

Vector Vector::operator-(const Vector& v2) {
    Vector v3;
    for (int i = 0; i < N1; ++i) {
        v3.v[i] = v[i] - v2.v[i];
    }
    return v3;
}

Vector& Vector::operator*(const Scalar& x) {
    for (int i = 0; i < N1; ++i) {
        v[i] *= x.get();
    }
    return *this;
}

double Vector::operator*(const Vector& v2) {
    double prod = 0;
    for (int i = 0; i < N1; ++i) {
        prod += v[i] * v2.v[i];
    }
    return prod;
}

Vector Vector::operator^(const Vector& v2) const {
    if (N1 != 3) throw std::invalid_argument("Cross product is only defined for 3D vectors");
    Vector v3;
    v3.v[0] = v[1] * v2.v[2] - v[2] * v2.v[1];
    v3.v[1] = v[2] * v2.v[0] - v[0] * v2.v[2];
    v3.v[2] = v[0] * v2.v[1] - v[1] * v2.v[0];
    return v3;
}

Vector& Vector::operator+=(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        v[i] += v2.v[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        v[i] -= v2.v[i];
    }
    return *this;
}

Vector& Vector::operator*=(const Scalar& x) {
    for (int i = 0; i < N1; ++i) {
        v[i] *= x.get();
    }
    return *this;
}

Vector& Vector::operator/=(const Scalar& x) {
    for (int i = 0; i < N1; ++i) {
        v[i] /= x.get();
    }
    return *this;
}

ntype& Vector::operator()(int index) {
    return v[index];
}

Vector& Vector::operator,(int val) {
    if (current_index < N1) {
        v[current_index++] = val;
    }
    return *this;
}

bool Vector::operator==(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        if (v[i] != v2.v[i]) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Vector& vec) {
    os << "[";
    for (int i = 0; i < N1; ++i) {
        os << vec.v[i];
        if (i != N1 - 1) os << ", ";
    }
    os << "]";
    return os;
}

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
    int d = 0;
    for (auto el : diag) {
        M[d][d] = el;
        d++;
    }
}

Matrix::~Matrix() {}

void Matrix::show(std::string name) const {
    std::cout << name << std::endl;
    for (int i = 0; i < N1; ++i) {
        std::cout << "|";
        for (int j = 0; j < N2; ++j) {
            std::cout << M[i][j];
            if (i < N1 - 1) 
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

Matrix Matrix::operator+(const Matrix& v2) {
    Matrix v3;
    for (int i = 0; i < N1; ++i) {
        v3.v[i] = v[i] + v2.v[i];
    }
    return v3;
}

Matrix Matrix::operator-(const Matrix& v2) {
    Matrix v3;
    for (int i = 0; i < N1; ++i) {
        v3.v[i] = v[i] - v2.v[i];
    }
    return v3;
}

ntype& Matrix::operator()(int i, int j) {
        return M[i][j];
    }