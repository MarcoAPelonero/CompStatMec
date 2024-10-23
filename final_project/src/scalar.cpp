#include "scalar.hpp"

Scalar::Scalar() : x(0) {}

Scalar::Scalar(ntype num) : x(num) {}

Scalar::~Scalar() {} 

ntype Scalar::get() const {
    return x;
}