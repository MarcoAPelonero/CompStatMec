#include "initialConfiguration.hpp"

initialConfiguration::initialConfiguration() {
    ensemble = particleEnsemble();
}

initialConfiguration::initialConfiguration(int N, ntype L) {
    ensemble = particleEnsemble(N, L);
}

initialConfiguration::~initialConfiguration() {}

void initialConfiguration::initializeSquareLattice(ntype L) {
    int numParticles = ensemble.getNumParticles();
    int latticeSpacing = rint(std::cbrt(numParticles));
    ntype spacing = L / latticeSpacing;

    while (latticeSpacing * latticeSpacing * latticeSpacing < numParticles) {
        latticeSpacing++;
        spacing = L / latticeSpacing;
    }

    int particleIndex = 0;
    for (int i = 0; i < latticeSpacing && particleIndex < numParticles; i++) {
        for (int j = 0; j < latticeSpacing && particleIndex < numParticles; j++) {
            for (int k = 0; k < latticeSpacing && particleIndex < numParticles; k++) {
                if (particleIndex < numParticles) {
                    ntype x = i * spacing;
                    ntype y = j * spacing;
                    ntype z = k * spacing;
                    ensemble.initParticle(particleIndex, Vector{x, y, z});
                    particleIndex++;
                }
            }
        }
    }
    ensemble.calculateEnergy();
}

particleEnsemble initialConfiguration::returnParticleEnsemble() {
    return ensemble;
}