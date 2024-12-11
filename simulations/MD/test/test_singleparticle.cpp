#include "particleEnsemble.hpp"

int main () {
    particleEnsemble ensemble(100, 10.0);
    ensemble.show();
    for (int i = 0; i < 20000; ++i) {
        ensemble.ensembleStep(0.1);
    }
    // ensemble.show();
    return 0;
}