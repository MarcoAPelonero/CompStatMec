#include "vector.hpp"

Vector::Vector() {
    for (int i = 0; i < N1; ++i) {
        v[i] = 0;
    }
}

Vector::Vector(ntype val) {
    for (int i = 0; i < N1; ++i) {
        v[i] = val;
    }
}

Vector::Vector(std::initializer_list<ntype> list) {
    int c = 0;
    for (auto el : list) {
        if (c < N1) {
            v[c++] = el;
        } else {
            break;
        }
    }
    for (; c < N1; ++c) {
        v[c] = 0.0;
    }
}

Vector::~Vector() {} 

void Vector::show(std::string name) const {
    std::cout << name << "(";
    for (int i = 0; i < N1; ++i) {
        std::cout << v[i];
        if (i < N1 - 1) 
            std::cout << ",";
    }
    std::cout << ")\n";
}

ntype Vector::get(int i) const {
    return v[i];
}

ntype Vector::set(int i, ntype val) {
    return v[i] = val;
}

double Vector::modulus() {
    return (*this) * (*this);
}

Vector random(ntype L) {
    std::random_device rd; 
    std::mt19937 eng(rd()); 

    std::uniform_real_distribution<double> distr(0.0, L);
    
    Vector v1;
    for (int i = 0; i < N1; i++) 
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
    for (int i = 0; i < N1; ++i) {
        v[i] = v2.v[i];
    }
    return *this;
}

Vector Vector::operator+(const Vector& v2) {
    Vector v3;
    for (int i = 0; i < N1; ++i) {
        v3.v[i] = v[i] + v2.v[i];
    }
    return v3;
}

Vector Vector::operator-(const Vector& v2) {
    Vector v3;
    for (int i = 0; i < N1; ++i) {
        v3.v[i] = v[i] - v2.v[i];
    }
    return v3;
}


Vector& Vector::operator*(const Scalar& x) {
    for (int i = 0; i < N1; ++i) {
        v[i] *= x.get();
    }
    return *this;
}

Vector operator*(const Scalar& x, const Vector& v1) {
    Vector v2;
    for (int i = 0; i<N1; i++)
        v2.v[i] = v1.v[i] * x.get();
    return v2;
} 

double Vector::operator*(const Vector& v2) {
    double prod = 0;
    for (int i = 0; i < N1; ++i) {
        prod += v[i] * v2.v[i];
    }
    return prod;
}

Vector Vector::operator^(const Vector& v2) const {
    if (N1 != 3) throw std::invalid_argument("Cross product is only defined for 3D vectors");
    Vector v3;
    v3.v[0] = v[1] * v2.v[2] - v[2] * v2.v[1];
    v3.v[1] = v[2] * v2.v[0] - v[0] * v2.v[2];
    v3.v[2] = v[0] * v2.v[1] - v[1] * v2.v[0];
    return v3;
}

Vector& Vector::operator+=(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        v[i] += v2.v[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        v[i] -= v2.v[i];
    }
    return *this;
}

Vector& Vector::operator*=(const Scalar& x) {
    for (int i = 0; i < N1; ++i) {
        v[i] *= x.get();
    }
    return *this;
} 

Vector& Vector::operator/=(const Scalar& x) {
    for (int i = 0; i < N1; ++i) {
        v[i] /= x.get();
    }
    return *this;
} 

ntype& Vector::operator()(int index) {
    return v[index];
}

Vector& Vector::operator,(int val) {
    if (current_index < N1) {
        v[current_index++] = val;
    }
    return *this;
}

bool Vector::operator==(const Vector& v2) {
    for (int i = 0; i < N1; ++i) {
        if (v[i] != v2.v[i]) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Vector& vec) {
    os << "[";
    for (int i = 0; i < N1; ++i) {
        os << vec.v[i];
        if (i != N1 - 1) os << ", ";
    }
    os << "]";
    return os;
}

Vector& Vector::mulcw(const Vector& v2) {
    for (int i = 0; i < N1; ++i) 
        this->v[i] = this->v[i] * v2.v[i];
    return *this;
}

Vector& Vector::divcw(const Vector& v2) {
    for (int i = 0; i < N1; ++i) 
        this->v[i] = this->v[i] / v2.v[i];
    return *this;
}

Vector mulcw(const Vector& v1, const Vector& v2) {
    Vector v3;
    for (int i = 0; i < N1; ++i)
        v3.v[i] = v1.v[i] * v2.v[i];
    return v3;
}

Vector divcw(const Vector& v1, const Vector& v2) {
    Vector v3;
    for (int i = 0; i < N1; ++i)
        v3.v[i] = v1.v[i] / v2.v[i];
    return v3;
}