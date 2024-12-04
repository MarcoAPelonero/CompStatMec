#include "initialConfiguration.hpp"

initialConfiguration::initialConfiguration() {
    ensemble = particleEnsemble();
}

initialConfiguration::initialConfiguration(int N){
    ensemble = particleEnsemble(N); 
}

initialConfiguration::~initialConfiguration() {}

void initialConfiguration::initializeSquareLattice(ntype L) {
    // Compute lattice spacing in each direction
    int numParticles = ensemble.getNumParticles();
    int latticeSpacing = rint(std::cbrt(numParticles));
    ntype spacing = L / latticeSpacing;

    // Adjust lattice spacing if the number of particles is not a perfect cube
    while (latticeSpacing * latticeSpacing * latticeSpacing < numParticles) {
        latticeSpacing++;
        spacing = L / latticeSpacing;
    }

    // Initialize particles in a square lattice
    int particleIndex = 0;
    for (int i = 0; i < latticeSpacing && particleIndex < numParticles; i++) {
        for (int j = 0; j < latticeSpacing && particleIndex < numParticles; j++) {
            for (int k = 0; k < latticeSpacing && particleIndex < numParticles; k++) {
                if (particleIndex < numParticles) {
                    // Set particle position based on lattice spacing
                    ntype x = i * spacing;
                    ntype y = j * spacing;
                    ntype z = k * spacing;
                    ensemble.initParticle(particleIndex, Vector{x,y,z});
                    particleIndex++;
                }
            }
        }
    }
}


particleEnsemble initialConfiguration::returnParticleEnsemble() {
    return ensemble;
}