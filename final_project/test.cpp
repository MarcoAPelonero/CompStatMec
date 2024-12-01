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
    
    // Default constructor
    Matrix mat1;
    mat1.show("mat1 (default constructor)");

    // Constructor with initializer list for diagonal
    Matrix mat2({1, 2, 3});
    mat2.show("mat2 (diagonal initializer)");

    // Constructor with rows and columns (multiplication of values)
    Matrix mat3({1, 2, 3}, {4, 5, 6});
    mat3.show("mat3 (rows and columns initializer)");

    // Set specific values
    mat1.set(0, 0, 5.0);
    mat1.set(1, 1, 10.0);
    std::cout << "Value at (0,0): " << mat1.get(0, 0) << std::endl;
    mat1.show("mat1 after using set");

    Matrix mat4({
    {4, 2, 0},
    {2, 4, 2},
    });
    mat4.show("mat4 before diagonalization");

    int swaps;
    Matrix diag_mat = mat4.diagonalize(swaps);
    diag_mat.show("mat4 after diagonalization");
    std::cout << "Number of swaps: " << swaps << std::endl;

    double det = mat4.determinant();
    std::cout << "Determinant of mat4: " << det << std::endl;

    // Assignment operator
    Matrix mat5;
    mat5 = mat4;
    mat5.show("mat5 (assigned from mat4)");

    // Addition operator
    Matrix mat6 = mat4 + mat5;
    mat6.show("mat6 (mat4 + mat5)");

    // Subtraction operator
    Matrix mat7 = mat4 - mat5;
    mat7.show("mat7 (mat4 - mat5)");

    // Multiplication operator (when you implement it)
    Matrix mat8 = mat4 * mat5;
    mat8.show("mat8 (mat4 * mat5)");

    // Operator() to access elements
    mat4(0, 0) = 100;
    std::cout << "Value at (0, 0) after using operator(): " << mat4(0, 0) << std::endl;
    mat4.show("mat4 after using operator()");

    std::cout << "New Vector class tests: " << mat4(0, 0) << std::endl;

    Vector vt = random_orient();
    vt.show("Random orient:");

    vt = random(10.0);


    return 0;
}
