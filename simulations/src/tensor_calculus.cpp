#include "tensor_calculus.hpp"

Vector operator+(Vector& v1, Scalar x) {
    Vector v2;
    for (int i = 0; i < N1; i++) 
        v2(i) = v1(i) + x.get();
    return v2;
}

Vector operator-(Vector& v1, Scalar x) {
    Vector v2;
    for (int i = 0; i < N1; i++) 
        v2(i) = v1(i) - x.get();
    return v2;
}

Vector operator*(Vector& v1, Scalar x) {
    Vector v2;
    for (int i = 0; i < N1; i++) 
        v2(i) = v1(i) * x.get();
    return v2;
}
Vector operator/(Vector& v1, Scalar x) {
    Vector v2;
    for (int i = 0; i < N1; i++) 
        v2(i) = v1(i) / x.get();
    return v2;
}
Vector& operator+=(Vector& v1, Scalar x) {
    for (int i = 0; i < N1; i++) 
        v1(i) += x.get();
}
Vector& operator-=(Vector& v1, Scalar x) {
    for (int i = 0; i < N1; i++) 
        v1(i) -= x.get();
}
Vector& operator*=(Vector& v1, Scalar x) {
    for (int i = 0; i < N1; i++) 
        v1(i) *= x.get();
}
Vector& operator/=(Vector& v1, Scalar x){
    for (int i = 0; i < N1; i++) 
        v1(i) /= x.get();
}

Vector operator*(Matrix& M1, Vector v1) {
    Vector v2;
    for (int i = 0; i < N1; i++) 
        for (int j = 0; j < N2; j++)
            v2(i) += M1(j,i) * v1(i);
    return v2;
}