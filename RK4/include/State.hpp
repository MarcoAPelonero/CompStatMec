#ifndef STATE_HPP
#define STATE_HPP

#include "Vec3.hpp"

struct State {
    Vec3 x;  // position
    Vec3 v;  // velocity
    Vec3 e;  // body axis unit vector
    Vec3 L;  // angular momentum in world frame
};

struct Deriv {
    Vec3 dx;
    Vec3 dv;
    Vec3 de;
    Vec3 dL;
};

#endif // STATE_HPP
