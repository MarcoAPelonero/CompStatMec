#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "vec.hpp"
#include <vector>
#include <fstream>
#include <iostream>

class Particle {
    private:
        Vector r;  // Position vector
        Vector rold; // Old position
        Vector v;
        Vector vold;  // Velocity vector
        Vector a;
        Vector aold; // Acceleration vector
        ntype m;  // Mass
    public:
        Particle();
        Particle(Vector r);
        Particle(Vector r, Vector v, Vector a, ntype m);
        Particle(const Particle& other); // Copy constructor
        ~Particle();

        void store();
        void restore();
        void random(ntype L);
        void setPosition(Vector pos);
        void setVelocity(Vector vel);
        void setOldPosition(Vector oldPos);\

        Vector getPosition();
        Vector getOldPosition();
        Vector getVelocity();
        Vector getOldVelocity();
        Vector getAcceleration();
        Vector getOldAcceleration();
        ntype getMass();
        
        Particle& operator=(Particle p);
};

#endif // PARTICLES_HPP