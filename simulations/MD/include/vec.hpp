#ifndef VEC_HPP
#define VEC_HPP

#include "randNumGen.hpp"
#include <initializer_list>
#include <string>
#include <iostream>
#include <random>

#define dim 3
#define ntype double

class Vector {
private:
    ntype r[dim];
public:
    Vector();
    Vector(ntype val);
    Vector(std::initializer_list<ntype> list);
    Vector(const Vector& other);
    ~Vector();
    ntype set(int i, ntype val);
    double modulus() const;
    Vector norm() const; 
    void random();
    Vector& random(ntype L);
    friend Vector random(ntype L);
    static Vector randomVector();

    void show(std::string name = "") const;
    ntype get(int i) const;

    Vector& operator=(const Vector& v2);
    Vector operator+(const Vector& v2) const;
    Vector operator-(const Vector& v2) const; // Added const version
    double operator*(const Vector& v2) const;
    Vector operator*(ntype scalar) const;
    Vector operator/(ntype scalar) const;

    Vector operator^(const Vector& v2) const;
    Vector& operator+=(const Vector& v2);
    Vector& operator-=(const Vector& v2);
    ntype& operator()(int index);
    const ntype& operator()(int index) const; // Added const version
    bool operator==(const Vector& v2) const;
    bool operator!=(const Vector& v2) const;
};

#endif // VEC_HPP