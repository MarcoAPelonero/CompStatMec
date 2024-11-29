#include "vec.hpp"

Vector::Vector() {
    for (int i = 0; i < dim; ++i) {
        r[i] = 0;
    }
}

Vector::Vector(ntype val) {
    for (int i = 0; i < dim; ++i) {
        r[i] = val;
    }
}

Vector::Vector(std::initializer_list<ntype> list) {
    int c = 0;
    for (auto el : list) {
        if (c < dim) {
            r[c++] = el;
        } else {
            break;
        }
    }
    for (; c < dim; ++c) {
        r[c] = 0.0;
    }
}

Vector::~Vector() {} 

void Vector::show(std::string name) const {
    std::cout << name << "(";
    for (int i = 0; i < dim; ++i) {
        std::cout << r[i];
        if (i < dim - 1) 
            std::cout << ",";
    }
    std::cout << ")\n";
}

ntype Vector::get(int i) const {
    return r[i];
}

ntype Vector::set(int i, ntype val) {
    return r[i] = val;
}

double Vector::modulus() {
    return (*this) * (*this);
}

Vector& Vector::random(ntype L) {
    std::random_device rd; 
    std::mt19937 eng(rd()); 

    std::uniform_real_distribution<double> distr(0.0, L);
    
    Vector v1;
    for (int i = 0; i < dim; i++) 
        v1(i) = distr(eng);
    return v1;
}

Vector random(ntype L) {
    std::random_device rd; 
    std::mt19937 eng(rd()); 

    std::uniform_real_distribution<double> distr(0.0, L);
    
    Vector v1;
    for (int i = 0; i < dim; i++) 
        v1(i) = distr(eng);
    return v1;
}

Vector random_orient() {
    std::random_device rd; 
    std::mt19937 eng(rd()); 

    std::uniform_real_distribution<double> distr(-1.0, 1.0);
    double a = 1, b = 1;

    while (a*a+b*b >= 1) {
        a = distr(eng);
        b = distr(eng);
    }
    
    double s = a*a+b*b;
    double x = 2 * a * sqrt(1 - s);
    double y = 2 * b * sqrt(1 - s);
    double z = 1 - 2 * s;

    return Vector {x,y,z};
} 

Vector& Vector::operator=(const Vector& v2) {
    for (int i = 0; i < dim; ++i) {
        r[i] = v2.r[i];
    }
    return *this;
}

Vector Vector::operator+(const Vector& v2) {
    Vector v3;
    for (int i = 0; i < dim; ++i) {
        v3.r[i] = r[i] + v2.r[i];
    }
    return v3;
}

Vector Vector::operator-(const Vector& v2) {
    Vector v3;
    for (int i = 0; i < dim; ++i) {
        v3.r[i] = r[i] - v2.r[i];
    }
    return v3;
}

double Vector::operator*(const Vector& v2) {
    double prod = 0;
    for (int i = 0; i < dim; ++i) {
        prod += r[i] * v2.r[i];
    }
    return prod;
}

Vector Vector::operator^(const Vector& v2) const {
    if (dim != 3) throw std::invalid_argument("Cross product is only defined for 3D vectors");
    Vector v3;
    v3.r[0] = r[1] * v2.r[2] - r[2] * v2.r[1];
    v3.r[1] = r[2] * v2.r[0] - r[0] * v2.r[2];
    v3.r[2] = r[0] * v2.r[1] - r[1] * v2.r[0];
    return v3;
}

Vector& Vector::operator+=(const Vector& v2) {
    for (int i = 0; i < dim; ++i) {
        r[i] += v2.r[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector& v2) {
    for (int i = 0; i < dim; ++i) {
        r[i] -= v2.r[i];
    }
    return *this;
}

ntype& Vector::operator()(int index) {
    return r[index];
}
