#ifndef VEC3_HPP
#define VEC3_HPP

#define _USE_MATH_DEFINES
#include <cmath>

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

inline double norm2(const Vec3& v) { return dot(v, v); }
inline double norm(const Vec3& v) { return std::sqrt(norm2(v)); }

inline Vec3 normalize(const Vec3& v, double eps = 1e-12) {
    double n = norm(v);
    if (n < eps) return Vec3(0,0,0);
    return v / n;
}

inline Vec3 clamp_unit(const Vec3& v, double eps = 1e-12) {
    Vec3 u = normalize(v, eps);
    if (norm(u) < eps) return Vec3(1,0,0);
    return u;
}

inline void split_parallel_perp(const Vec3& v, const Vec3& e_unit, Vec3& v_par, Vec3& v_perp) {
    v_par = e_unit * dot(v, e_unit);
    v_perp = v - v_par;
}

inline double sin_alpha(const Vec3& v, const Vec3& e_unit, double eps = 1e-12) {
    double vmag = norm(v);
    if (vmag < eps) return 0.0;
    Vec3 vhat = v / vmag;
    double c = dot(vhat, e_unit);
    double s2 = 1.0 - c*c;
    return (s2 > 0.0) ? std::sqrt(s2) : 0.0;
}

#endif // VEC3_HPP
