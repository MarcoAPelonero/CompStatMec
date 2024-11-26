#ifndef SIM_HPP
#define SIM_HPP

#include "config.hpp"
#include "particle.hpp"
#include <fstream>
#include <cmath>
#include "function.hpp"

class sys {
    private:
        std::vector<Particle> old_particles;
        std::vector<Particle> new_particles;
    public:
        sys() {
            old_particles = std::vector<Particle>(N_PARTICLES);
            new_particles = std::vector<Particle>(N_PARTICLES);
        }

        void store() {
            for (int i = 0; i < N_PARTICLES; ++i) {
                old_particles[i] = new_particles[i];
                new_particles[i].store();
            }
        }

        void restore() {
            for (int i = 0; i < N_PARTICLES; ++i) {
                new_particles[i] = old_particles[i];
                new_particles[i].restore();
            }
        }

        void show() {
            for (int i = 0; i < N_PARTICLES; ++i) {
                std::cout << "Particle " << i << ":\n";
                std::cout << "Position: ";
                new_particles[i].getPosition().show();
                std::cout << "Sigma: " << new_particles[i].getSigma() << "\n";
                std::cout << "Epsilon: " << new_particles[i].getEpsilon() << "\n";
                std::cout << "Rcut: " << new_particles[i].getRcut() << "\n";
            }
        }
};

class simulation: public sys {
    
    public:
        void prepare_configuration();
        void update_configuration();
        void init_rng(int n);   
        void run(void) {
            // intentionally void
        };
};

#endif