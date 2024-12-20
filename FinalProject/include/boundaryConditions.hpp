#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "particle.hpp"

class BoundaryConditions {
private:
    double boxLength;

public:
    BoundaryConditions(double L);
    void apply(Particle &p) const;
};

#endif // BOUNDARY_CONDITIONS_HPP
