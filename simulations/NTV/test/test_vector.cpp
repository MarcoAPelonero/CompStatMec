#include <iostream>
#include "vec.hpp"

void testDefaultConstructor() {
    Vector v;
    std::cout << "Default constructor test passed." << std::endl;
}

void testValueConstructor() {
    Vector v(5.0);
    for (int i = 0; i < dim; ++i) {
        if (v.get(i) != 5.0) {
            std::cerr << "Value constructor test failed." << std::endl;
            return;
        }
    }
    std::cout << "Value constructor test passed." << std::endl;
}

void testInitializerListConstructor() {
    Vector v{1.0, 2.0, 3.0};
    if (v.get(0) == 1.0 && v.get(1) == 2.0 && v.get(2) == 3.0) {
        std::cout << "Initializer list constructor test passed." << std::endl;
    } else {
        std::cerr << "Initializer list constructor test failed." << std::endl;
    }
}

void testAssignmentOperator() {
    Vector v1{1.0, 2.0, 3.0};
    Vector v2 = v1;
    if (v2.get(0) == 1.0 && v2.get(1) == 2.0 && v2.get(2) == 3.0) {
        std::cout << "Assignment operator test passed." << std::endl;
    } else {
        std::cerr << "Assignment operator test failed." << std::endl;
    }
}

void testAdditionOperator() {
    Vector v1{1.0, 2.0, 3.0};
    Vector v2{4.0, 5.0, 6.0};
    Vector v3 = v1 + v2;
    if (v3.get(0) == 5.0 && v3.get(1) == 7.0 && v3.get(2) == 9.0) {
        std::cout << "Addition operator test passed." << std::endl;
    } else {
        std::cerr << "Addition operator test failed." << std::endl;
    }
}

void testSubtractionOperator() {
    Vector v1{4.0, 5.0, 6.0};
    Vector v2{1.0, 2.0, 3.0};
    Vector v3 = v1 - v2;
    if (v3.get(0) == 3.0 && v3.get(1) == 3.0 && v3.get(2) == 3.0) {
        std::cout << "Subtraction operator test passed." << std::endl;
    } else {
        std::cerr << "Subtraction operator test failed." << std::endl;
    }
}

void testDotProductOperator() {
    Vector v1{1.0, 2.0, 3.0};
    Vector v2{4.0, 5.0, 6.0};
    double dotProduct = v1 * v2;
    if (dotProduct == 32.0) {
        std::cout << "Dot product operator test passed." << std::endl;
    } else {
        std::cerr << "Dot product operator test failed." << std::endl;
    }
}

void testCrossProductOperator() {
    Vector v1{1.0, 2.0, 3.0};
    Vector v2{4.0, 5.0, 6.0};
    Vector v3 = v1 ^ v2;
    if (v3.get(0) == -3.0 && v3.get(1) == 6.0 && v3.get(2) == -3.0) {
        std::cout << "Cross product operator test passed." << std::endl;
    } else {
        std::cerr << "Cross product operator test failed." << std::endl;
    }
}

void testScalarMultiplicationOperator() {
    Vector v1{1.0, 2.0, 3.0};
    Vector v2 = v1 * 2.0;
    if (v2.get(0) == 2.0 && v2.get(1) == 4.0 && v2.get(2) == 6.0) {
        std::cout << "Scalar multiplication operator test passed." << std::endl;
    } else {
        std::cerr << "Scalar multiplication operator test failed." << std::endl;
    }
}

void testModulus() {
    Vector v{3.0, 4.0, 0.0};
    double mod = v.modulus();
    if (mod == 5.0) {
        std::cout << "Modulus test passed." << std::endl;
    } else {
        std::cerr << "Modulus test failed." << std::endl;
    }
}

void testRandom() {
    Vector v;
    v.random();
    bool passed = true;
    for (int i = 0; i < dim; ++i) {
        if (v.get(i) < 0.0 || v.get(i) > 1.0) {
            passed = false;
            break;
        }
    }
    if (passed) {
        std::cout << "Random test passed." << std::endl;
    } else {
        std::cerr << "Random test failed." << std::endl;
    }
}

void testRandomWithL() {
    Vector v;
    v.random(10.0);
    bool passed = true;
    for (int i = 0; i < dim; ++i) {
        if (v.get(i) < 0.0 || v.get(i) > 10.0) {
            passed = false;
            break;
        }
    }
    if (passed) {
        std::cout << "Random with L test passed." << std::endl;
    } else {
        std::cerr << "Random with L test failed." << std::endl;
    }
}

void testShow() {
    Vector v{1.0, 2.0, 3.0};
    v.show("Vector v");
}

int main() {
    testDefaultConstructor();
    testValueConstructor();
    testInitializerListConstructor();
    testAssignmentOperator();
    testAdditionOperator();
    testSubtractionOperator();
    testDotProductOperator();
    testCrossProductOperator();
    testScalarMultiplicationOperator();
    testModulus();
    testRandom();
    testRandomWithL();
    testShow();
    return 0;
}