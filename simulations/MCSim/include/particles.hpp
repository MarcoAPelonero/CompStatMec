#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "vec.hpp"
#include "configuration.hpp"
#include "interactionPotentials.hpp"

// Forward declarations to avoid circular dependencies
class initialConfiguration;

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
        ~Particle();

        void store();
        void restore();
        void random(ntype L);
        void setPosition(Vector pos);
        void setVelocity(Vector vel);
        Vector getPosition();
        Vector getOldPosition();

        Particle& operator=(Particle p);
};

class particleEnsemble: public interactionPotential {  // Removed initialConfiguration from inheritance
    private:
        std::vector<Particle> particles;
        int numParticles;  // Number of particles
        ntype ensembleEnergy;
        initialConfiguration* config;  // Member instead of inheritance
    public:
        particleEnsemble();
        particleEnsemble(int N);
        ~particleEnsemble();

        Particle& operator()(int i);
        int getNumParticles();

        void initializeRandom(ntype L);
        void initParticle(int i, Vector v);
        void store();
        void restore();

        ntype getEnergy();
        ntype calculateEnergy();
        void updateEnergy(ntype newEnergy);
};

#endif // PARTICLE_HPP
