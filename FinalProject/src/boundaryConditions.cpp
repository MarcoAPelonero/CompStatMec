#include "boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(double L)
    : boxLength(L)
{}

void BoundaryConditions::apply(Particle &p) const {
    Vec pos = p.getPosition();
    double x = pos.getX();
    double y = pos.getY();
    double z = pos.getZ();

    if(x < 0.0) x += boxLength;
    if(x >= boxLength) x -= boxLength;
    if(y < 0.0) y += boxLength;
    if(y >= boxLength) y -= boxLength;
    if(z < 0.0) z += boxLength;
    if(z >= boxLength) z -= boxLength;

    p.setPosition(Vec(x, y, z));
}
