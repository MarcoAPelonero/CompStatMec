#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "vec.hpp"
#include "configuration.hpp"
#include "interactionPotentials.hpp"
#include <vector>
#include <fstream>
#include <iostream>

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
        void setOldPosition(Vector oldPos);
        Vector getPosition();
        Vector getOldPosition();
        Vector getVelocity();
        Vector getOldVelocity();

        Particle& operator=(Particle p);
};

class particleEnsemble: public interactionPotential {  // Removed initialConfiguration from inheritance
    private:
        std::vector<Particle> particles;
        int numParticles;  // Number of particles
        ntype ensembleEnergy;
        initialConfiguration* config;  // Member instead of inheritance
        ntype boxLength;
    public:
        particleEnsemble();
        particleEnsemble(int N,ntype L);
        ~particleEnsemble();

        Particle& operator()(int i);
        int getNumParticles();

        void initializeRandom(ntype L);
        void initParticle(int i, Vector v);
        void store();
        void restore();

        ntype getEnergy();
        ntype calculateEnergy();
        void setInitialEnergy();
        // void particleEnsemble::setEnergy(ntype newEnergy);
        void updateEnergy(ntype deltaEnergy);

        void saveSnapshotPosition(std::ofstream& file, int step);
        void saveSnapshotVelocity(std::ofstream& file);
        void ensembleSnapshot(std::ofstream& file); 
};

#endif // PARTICLE_HPP