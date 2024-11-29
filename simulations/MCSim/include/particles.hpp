#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "vec.hpp"
#include "config.hpp"
#include <vector>

class Particle {
    private:
        Vector r;  // Position vector
        Vector rold; //Old position
        Vector v;
        Vector vold;  // Velocity vector
        Vector a;
        Vector aold; // Acceleration vector
        ntype m;  // Mass
    public:
        Particle();
        Particle(Vector r);
        Particle(Vector r, Vector v, Vector a, ntype m);
        ~Particle();

        void store();
        void restore();
        void random();
        Vector getPosition();
};


class particleEnsemble: public Particle {
    private:
        std::vector<Particle> particles;
    public:
        particleEnsemble(int N);
        ~particleEnsemble();

        Particle& operator()(int i);

        void initialize(ntype L);
        void initializeRandom(ntype L);
        void store();
        void restore();
        void tra_move(Vector delr);
};
#endif // PARTICLE_HPP 