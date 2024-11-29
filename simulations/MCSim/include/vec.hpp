#ifndef VEC_HPP
#define VEC_HPP

#include "config.hpp"
#include "randNumGen.hpp"
#include <initializer_list>
#include <string>
#include <iostream>
#include <random>

class Vector {
    private:
        ntype r[dim];
        rng.seed(42);
    public:
    Vector();
    Vector(ntype val);
    Vector(std::initializer_list<ntype> list);
    ~Vector(); 

    void show(std::string name = "") const;
    ntype get(int i) const;
    ntype set(int i, ntype val);
    double modulus();
    Vector vretint();
    Vector& random(ntype L);
    friend Vector random(ntype L);
    friend Vector random_orient();

    Vector& operator=(const Vector& v2);
    Vector operator+(const Vector& v2);
    Vector operator-(const Vector& v2);
    
    double operator*(const Vector& v2);
    Vector operator^(const Vector& v2) const; 
    Vector& operator+=(const Vector& v2);
    Vector& operator-=(const Vector& v2);
    ntype& operator()(int index);
    Vector& operator,(int val);
    bool operator==(const Vector& v2);
};

#endif // VEC_HPP