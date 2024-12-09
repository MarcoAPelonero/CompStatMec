#ifndef VEC_HPP
#define VEC_HPP

#include "randNumGen.hpp"
#include "configuration.hpp"
#include <initializer_list>
#include <string>
#include <iostream>
#include <random>

class Vector {
    private:
        double r[dim];
    public:
    Vector();
    Vector(ntype val);
    Vector(std::initializer_list<ntype> list);
    ~Vector(); 

    void show(std::string name = "") const;
    ntype get(int i) const;
    ntype set(int i, ntype val);
    double modulus();
    void random();
    Vector& random(ntype L);
    friend Vector random(ntype L);
    static Vector randomVector();

    Vector& operator=(const Vector& v2);
    Vector operator+(const Vector& v2);
    Vector operator-(const Vector& v2);
    
    double operator*(const Vector& v2);
    Vector operator*(ntype scalar);
    Vector operator^(const Vector& v2) const; 
    Vector& operator+=(const Vector& v2);
    Vector& operator-=(const Vector& v2);
    ntype& operator()(int index);
};

#endif // VEC_HPP