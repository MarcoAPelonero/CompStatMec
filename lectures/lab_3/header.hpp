#include <cstdlib>
#include <iostream>

double a;
// Vediamo come si usa un template e perchè è molto fico
// è un parametro, dipendentemente a come assegno ntype, ho una certa classe

template <typename ntype>

class genvec {
    private:
        ntype v[3];
    public:
        void set(int i, ntype value) { v[i] = value;}
};
