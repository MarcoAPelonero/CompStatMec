#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "particleEnsemble.hpp"
#include "monteCarloMove.hpp"

class MonteCarloSimulation {
private:
    particleEnsemble ensemble;
    MonteCarloMove mcMove;
    double temperature;
    double displacement;
    int steps;  
public:
    MonteCarloSimulation(int numParticles, double boxLength, double T, double delta, int numSteps);
    ~MonteCarloSimulation();

    void initRNG();
    
    // Modified runSimulation to accept runId
    void runSimulation(int snapshotInterval, int runId);
    void printAcceptanceRatio() const;
    double getAcceptanceRatio() const;
    double getEnergy() const;   
};

#endif // MONTECARLOSIMULATION_HPP
