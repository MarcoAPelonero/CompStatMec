#ifndef PARTICLEENSEMBLE_HPP
#define PARTICLEENSEMBLE_HPP

#include "particle.hpp"
#include "interactionPotential.hpp"
#include "progressBar.hpp"
#include <vector>
#include <fstream>

class particleEnsemble {
    private:
        std::vector<Particle> particles;
        int numParticles;
        ntype boxSize;
        ntype ensembleEnergy;
        interactionPotential potential;

    public:
        particleEnsemble(int numParticles, ntype boxSize);
        ~particleEnsemble();

        Particle& operator()(int i);
        int getNumParticles();

        void show();
        void storeEnsemble();
        
        ntype getEnsembleEnergy();

        void applyPeriodicBoundary(Vector &pos);

        void ensembleSnapshot(std::ofstream &outFile, bool writeHeader = false);

        void particleEnergyStep(int particleIndex);
        void ensembleStep();
        void ensembleAverageEnergy();

        void Euler(ntype dt, int numSteps, std::ofstream &outFile);
        void EulerCromer(ntype dt, int numSteps, std::ofstream &outFile);
        void SpeedVerlet(ntype dt, int numSteps, std::ofstream &outFile); 
        void PredictorCorrector(ntype dt, int numSteps, std::ofstream &outFile);
        void ThermoSpeedVerlet(ntype dt, int numSteps, std::ofstream &outFile, ntype temperature = 1.0, int thermoSteps = 10);  
};

#endif // PARTICLEENSEMBLE_HPP
