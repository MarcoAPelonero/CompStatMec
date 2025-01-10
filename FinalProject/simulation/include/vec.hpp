#ifndef VEC_HPP
#define VEC_HPP

#include <iostream>
#include <cmath>

#define DIMENSION 3

class Vec {
    private:
        double x, y, z;
    public:
        Vec() : x(0), y(0), z(0) {}
        Vec(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
        Vec(const Vec &v) : x(v.x), y(v.y), z(v.z) {}

        double getX() const { return x; }
        double getY() const { return y; }
        double getZ() const { return z; }

        void setX(double X) { this->x = X; }
        void setY(double Y) { this->y = Y; }
        void setZ(double Z) { this->z = Z; }

        void show(const std::string &prefix = "") const {
            std::cout << prefix << "(" << x << ", " << y << ", " << z << ")" << std::endl;
        }

        double length() const {
            return sqrt(x * x + y * y + z * z);
        }

        Vec normalize() const {
            double l = length();
            return Vec(x / l, y / l, z / l);
        }

        double dot(const Vec &v) const {
            return x * v.x + y * v.y + z * v.z;
        }

        Vec cross(const Vec &v) const {
            return Vec(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
        }

        Vec& operator=(const Vec &v) {
            x = v.x;
            y = v.y;
            z = v.z;
            return *this;
        }

        Vec operator+(const Vec &v) const {
            return Vec(x + v.x, y + v.y, z + v.z);
        }
        Vec operator-(const Vec &v) const {
            return Vec(x - v.x, y - v.y, z - v.z);
        }
        Vec operator*(double s) const {
            return Vec(x * s, y * s, z * s);
        }
        Vec operator/(double s) const {
            return Vec(x / s, y / s, z / s);
        }
        Vec operator-() const {
            return Vec(-x, -y, -z);
        }

        Vec operator+=(const Vec &v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        Vec operator-=(const Vec &v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        Vec operator*=(double s) {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }

        Vec operator/=(double s) {
            x /= s;
            y /= s;
            z /= s;
            return *this;
        }

        bool operator==(const Vec &v) const {
            return x == v.x && y == v.y && z == v.z;
        }

        double operator()(int i) const {
            if (i == 0) return x;
            if (i == 1) return y;
            if (i == 2) return z;
            throw std::out_of_range("Index out of range in Vec::operator()");
        }

        friend Vec operator*(double s, const Vec &v) {
        return Vec(v.x * s, v.y * s, v.z * s);
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec& vec) {

        os << vec.getX() << " " << vec.getY() << " " << vec.getZ();

        return os;

    }
};

#endif // VEC_HPP