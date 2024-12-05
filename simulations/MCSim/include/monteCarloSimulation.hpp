#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "particles.hpp"
#include "monteCarloMove.hpp"
#include "interactionPotentials.hpp"
#include <fstream>
#include "monteCarloSimulation.hpp"
#include "progressBar.hpp"
#include <iostream>

class MonteCarloSimulation: public interactionPotential {
private:
    particleEnsemble ensemble;
    interactionPotential potential;
    MonteCarloMove mcMove;
    int numSteps;
    ntype L;

    ntype Energy;

    std::ofstream energyFile;
    std::ofstream positionFile;
    std::ofstream velocityFile;
public:
    MonteCarloSimulation(int N, ntype T, ntype L, int steps);
    ~MonteCarloSimulation();

    void initRNG();
    void run();
};

#endif // MONTECARLOSIMULATION_HPP