// C++
#include "particleEnsemble.hpp"
#include <cmath>
#include <random>

particleEnsemble::particleEnsemble(int numParticles, ntype boxSize) {
    this->numParticles = numParticles;
    this->boxSize = boxSize;

    int particlesPerSide = std::ceil(std::cbrt(numParticles));
    ntype spacing = boxSize / particlesPerSide;

    std::default_random_engine generator;
    std::normal_distribution<ntype> distribution(0.0, 1.0);

    int count = 0;
    for (int i = 0; i < particlesPerSide && count < numParticles; ++i) {
        for (int j = 0; j < particlesPerSide && count < numParticles; ++j) {
            for (int k = 0; k < particlesPerSide && count < numParticles; ++k) {
                Particle p(Vector{i * spacing, j * spacing, k * spacing});

                ntype vx = distribution(generator);
                ntype vy = distribution(generator);
                ntype vz = distribution(generator);
                p.setVelocity(Vector{vx, vy, vz});

                particles.push_back(p);
                ++count;
                if (count >= numParticles) break;
            }
        }
    }
}

void particleEnsemble::show() {
    for (int i = 0; i < numParticles; ++i) {
        std::cout << "Particle " << i << ":" << std::endl;
        particles[i].show();
    }
}

