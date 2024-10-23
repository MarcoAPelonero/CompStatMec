#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "config.hpp"
#include <initializer_list>
#include <string>
#include <iostream>

class Matrix {
    private:
        ntype M[N1][N2];
    public:
        Matrix();
        Matrix(ntype val);
        Matrix(std::initializer_list<ntype> rows, std::initializer_list<ntype> cols);
        Matrix(std::initializer_list<ntype> diag);
        ~Matrix(); 

        void show(std::string name = "") const;
        ntype get(int i, int j) const;
        ntype set(int i, int j, ntype val);
        Matrix diagonalize(int &swaps) const;
        double determinant() const;

        Matrix& operator=(const Matrix& M2);
        Matrix operator+(const Matrix& M2);
        Matrix operator-(const Matrix& M2);
        Matrix operator*(const Matrix& M2);
        ntype& operator()(int i, int j);
}; 

#endif