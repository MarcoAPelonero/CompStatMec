#include <iostream>
#include "header.hpp"
// Impariamo a usare i namespace

namespace nsA {
    int a;
    int b;
}

namespace nsB {
    int a;
    int b;
}

namespace nsC {
    int a;
    int b;
}

// Quando uso il mio namespace personale faccio un override di std, e faccio casini
//using namespace nsC;

int main() {
    // Questo è uno scoping 
    nsA::a = 1;
    nsB::a = 2;
    nsC::a = 3;
    std::cout << "A " << nsA::a << std::endl;
    std::cout << "C " << nsC::a << std::endl;
    //std::cout << "Cnew " << a << std::endl;
    return 0;
    // La cosa molto figa di namespace, è che non posso usare solo gli standard type, ma anche i tipi custom delle mie classi
    genvec<double> vec;
}