#ifndef TENSORS_HPP
#define TENSORS_HPP

#include "config.hpp"
#include <initializer_list>
#include <string>

class Scalar {
    private: 
        ntype x;
    public: 
        Scalar();
        Scalar(ntype num);

        ntype get() const;
};

class Vector {
    private: 
        ntype v[N1];
    public: 
        Vector();
        Vector(std::initializer_list<ntype> list);
        
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
};

#endif TENSORS_HPP