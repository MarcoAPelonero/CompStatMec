#ifndef TENSOR_CALCULUS_HPP
#define TENSOR_CALCULUS_HPP

#include "scalar.hpp"
#include "vector.hpp"
#include "matrix.hpp"

Vector operator+(Vector& v1, Scalar x);
Vector operator-(Vector& v1, Scalar x);
Vector operator*(Vector& v1, Scalar x);
Vector operator/(Vector& v1, Scalar x);
Vector& operator+=(Vector& v1, Scalar x);
Vector& operator-=(Vector& v1, Scalar x);
Vector& operator*=(Vector& v1, Scalar x);
Vector& operator/=(Vector& v1, Scalar x);

Vector operator*(Matrix& M1, Vector v1);

#endif