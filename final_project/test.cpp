// test.cpp
#include <iostream>
#include "Scalar.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

int main() {
    // Testing Scalar class
    Scalar s1;
    Scalar s2(5.0);
    std::cout << "Scalar s1 value: " << s1.get() << std::endl;
    std::cout << "Scalar s2 value: " << s2.get() << std::endl;

    // Testing Vector class
    Vector v1;
    v1.show("v1");

    Vector v2 = {1.0, 2.0, 3.0};
    v2.show("v2");

    // Testing get and set methods
    v1.set(0, 4.0);
    v1.set(1, 5.0);
    v1.set(2, 6.0);
    v1.show("v1 after set");

    std::cout << "v1[0]: " << v1.get(0) << std::endl;
    std::cout << "v1[1]: " << v1.get(1) << std::endl;
    std::cout << "v1[2]: " << v1.get(2) << std::endl;

    // Testing operator overloading
    Vector v3 = v1 + v2;
    v3.show("v3 = v1 + v2");

    Vector v4 = v1 - v2;
    v4.show("v4 = v1 - v2");

    Scalar s3(2.0);
    Vector v5 = v1;
    // v5 = v5 * s3;
    v5.show("v5 = v1 * s3");

    double dot_product = v1 * v2;
    std::cout << "Dot product of v1 and v2: " << dot_product << std::endl;

    Vector v6 = v1 ^ v2;
    v6.show("v6 = v1 ^ v2 (cross product)");

    v1 += v2;
    v1.show("v1 after v1 += v2");

    v1 -= v2;
    v1.show("v1 after v1 -= v2");
    /*
    v1 *= s3;
    v1.show("v1 after v1 *= s3");

    v1 /= s3;
    v1.show("v1 after v1 /= s3");
    */
    // Testing modulus
    double modulus_v1 = v1.modulus();
    std::cout << "Modulus of v1: " << modulus_v1 << std::endl;

    // Testing operator()
    v1(0) = 10.0;
    v1(1) = 20.0;
    v1(2) = 30.0;
    v1.show("v1 after using operator()");

    // Testing equality operator
    bool are_equal = (v1 == v1);
    std::cout << "v1 == v1: " << are_equal << std::endl;

    are_equal = (v1 == v2);
    std::cout << "v1 == v2: " << are_equal << std::endl;

    return 0;
}
