#include <stdio.h>

class myvec {
    private:
        int *v;
        int len;
    public:
        // Quando creo l'ogetto e vuoto
        myvec () {
            len = 0;
            v=nullptr;
        }

        void size(int l) {
            len = l;
            // New e la funzione di cpp che fa allocazione dinamica, anche di un tipo custom
            v = new int[len];
        }
        
        void set (int i, int val) {
            v[i] = val;
        }

        int get (int i) const {
            return v[i];
        }

        // Altri metodi
};

