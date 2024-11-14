#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "config.hpp"
#include "scalar.hpp"
#include <initializer_list>
#include <string>
#include <iostream>
#include <random>

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
    friend Vector random(ntype L);
    friend Vector random_orient();

    Vector& operator=(const Vector& v2);
    Vector operator+(const Vector& v2);
    Vector operator-(const Vector& v2);
    friend Vector operator*(const Scalar& x, const Vector& v1);
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
    friend Vector mulcw(const Vector& v1, const Vector& v2);
    friend Vector divcw(const Vector& v1, const Vector& v2); 
    Vector& mulcw(const Vector& v2);
    Vector& divcw(const Vector& v2);

};

#endif