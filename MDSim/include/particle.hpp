#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "vec.hpp"

class Particle {
    private:
        double mass;
        Vec position, velocity, force;
        double kineticEnergy, potentialEnergy, totalEnergy;
    public:
        Particle() : mass(1.0), position(Vec()), velocity(Vec()), force(Vec()), kineticEnergy(0), potentialEnergy(0), totalEnergy(0) {}
        Particle(double m, const Vec &p, const Vec &v) : mass(m), position(p), velocity(v), force(Vec()), kineticEnergy(0), potentialEnergy(0), totalEnergy(0) {}

        double getMass() const { return mass; }
        Vec getPosition() const { return position; }
        Vec getVelocity() const { return velocity; }
        Vec getForce() const { return force; }
        double getKineticEnergy() const { return kineticEnergy; }
        double getPotentialEnergy() const { return potentialEnergy; }
        double getTotalEnergy() const { return totalEnergy; }

        void setMass(double m) { mass = m; }
        void setPosition(const Vec &p) { position = p; }
        void setVelocity(const Vec &v) { velocity = v; }
        void setForce(const Vec &f) { force = f; }
        void setKineticEnergy(double ke) { kineticEnergy = ke; }
        void setPotentialEnergy(double pe) { potentialEnergy = pe; }
        void setTotalEnergy(double te) { totalEnergy = te; }

        void addForce(const Vec &f) { force = force + f; }

        void show(const std::string &prefix = "") const {
            std::cout << prefix << "Mass: " << mass << std::endl;
            position.show(prefix + "Position: ");
            velocity.show(prefix + "Velocity: ");
            force.show(prefix + "Force: ");
            std::cout << prefix << "Kinetic energy: " << kineticEnergy << std::endl;
            std::cout << prefix << "Potential energy: " << potentialEnergy << std::endl;
            std::cout << prefix << "Total energy: " << totalEnergy << std::endl;
        }

        void computeKineticEnergy() {
            kineticEnergy = 0.5 * mass * velocity.length() * velocity.length();
        }

        void computeTotalEnergy() {
            totalEnergy = kineticEnergy + potentialEnergy;
        }
};

#endif // PARTICLE_HPP