#ifndef TENSORS_HPP
#define TENSORS_HPP

#include "config.hpp"
#include <initializer_list>
#include <string>
#include <iostream>

class Scalar {
private:
    ntype x;
public:
    Scalar();
    Scalar(ntype num);
    ~Scalar(); 

    ntype get() const;
};

class Vector {
private:
    ntype v[N1];
    int current_index = 0; 
public:
    Vector();
    Vector(ntype val);
    Vector(std::initializer_list<ntype> list);
    ~Vector(); 

    void show(std::string name = "") const;
    ntype get(int i) const;
    ntype set(int i, ntype val);
    double modulus();

    Vector& operator=(const Vector& v2);
    Vector operator+(const Vector& v2);
    Vector operator-(const Vector& v2);
    Vector& operator*(const Scalar& x);
    double operator*(const Vector& v2);
    Vector operator^(const Vector& v2) const; 
    Vector& operator+=(const Vector& v2);
    Vector& operator-=(const Vector& v2);
    Vector& operator*=(const Scalar& x);
    Vector& operator/=(const Scalar& x);
    ntype& operator()(int index);
    Vector& operator,(int val);
    bool operator==(const Vector& v2);

    friend std::ostream& operator<<(std::ostream& os, const Vector& v); // Declare friend functions
};

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

        ntype& operator()(int i, int j);
};      

#endif // TENSORS_HPP
