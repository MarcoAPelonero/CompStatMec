#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "config.hpp"
#include <initializer_list>
#include <string>
#include <iostream>

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
    // Vector& operator*(const Scalar& x);
    double operator*(const Vector& v2);
    Vector operator^(const Vector& v2) const; 
    Vector& operator+=(const Vector& v2);
    Vector& operator-=(const Vector& v2);
    // Vector& operator*=(const Scalar& x);
    // Vector& operator/=(const Scalar& x);
    ntype& operator()(int index);
    Vector& operator,(int val);
    bool operator==(const Vector& v2);

    friend std::ostream& operator<<(std::ostream& os, const Vector& v); // Declare friend functions
};

#endif