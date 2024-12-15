#include <iostream>
#include <cmath>
#include "vec.hpp"

void testDefaultConstructor() {
    Vec v;
    if (v.getX() == 0 && v.getY() == 0 && v.getZ() == 0) {
        std::cout << "Default constructor test passed." << std::endl;
    } else {
        std::cout << "Default constructor test failed." << std::endl;
    }
}

void testParameterizedConstructor() {
    Vec v(1.0, 2.0, 3.0);
    if (v.getX() == 1.0 && v.getY() == 2.0 && v.getZ() == 3.0) {
        std::cout << "Parameterized constructor test passed." << std::endl;
    } else {
        std::cout << "Parameterized constructor test failed." << std::endl;
    }
}

void testLength() {
    Vec v(3.0, 4.0, 0.0);
    if (std::abs(v.length() - 5.0) < 1e-9) {
        std::cout << "Length test passed." << std::endl;
    } else {
        std::cout << "Length test failed." << std::endl;
    }
}

void testNormalize() {
    Vec v(3.0, 4.0, 0.0);
    Vec n = v.normalize();
    if (std::abs(n.getX() - 0.6) < 1e-9 && std::abs(n.getY() - 0.8) < 1e-9 && std::abs(n.getZ()) < 1e-9) {
        std::cout << "Normalize test passed." << std::endl;
    } else {
        std::cout << "Normalize test failed." << std::endl;
    }
}

void testDotProduct() {
    Vec v1(1.0, 2.0, 3.0);
    Vec v2(4.0, -5.0, 6.0);
    if (std::abs(v1.dot(v2) - 12.0) < 1e-9) {
        std::cout << "Dot product test passed." << std::endl;
    } else {
        std::cout << "Dot product test failed." << std::endl;
    }
}

void testCrossProduct() {
    Vec v1(1.0, 2.0, 3.0);
    Vec v2(4.0, 5.0, 6.0);
    Vec v3 = v1.cross(v2);
    if (std::abs(v3.getX() + 3.0) < 1e-9 && std::abs(v3.getY() - 6.0) < 1e-9 && std::abs(v3.getZ() + 3.0) < 1e-9) {
        std::cout << "Cross product test passed." << std::endl;
    } else {
        std::cout << "Cross product test failed." << std::endl;
    }
}

void testShow() {
    Vec v(1.0, 2.0, 3.0);
    std::cout << "Expected output: Vector: (1, 2, 3)" << std::endl;
    std::cout << "Actual output: ";
    v.show("Vector: ");
}

void testFriendOperator() {
    Vec v(1.0, 2.0, 3.0);
    Vec result = 2.0 * v;
    if (result.getX() == 2.0 && result.getY() == 4.0 && result.getZ() == 6.0) {
        std::cout << "Friend operator test passed." << std::endl;
    } else {
        std::cout << "Friend operator test failed." << std::endl;
    }
}

void testNegationOperator() {
        Vec v(1.0, 2.0, 3.0);
        Vec result = -v;
        if (result.getX() == -1.0 && result.getY() == -2.0 && result.getZ() == -3.0) {
            std::cout << "Negation operator test passed." << std::endl;
        } else {
            std::cout << "Negation operator test failed." << std::endl;
        }
}

int main() {
    testDefaultConstructor();
    testParameterizedConstructor();
    testLength();
    testNormalize();
    testDotProduct();
    testCrossProduct();
    testShow();
    testFriendOperator();
    testNegationOperator();
    return 0;
}