#ifndef RANDNUMGEN_HPP
#define RANDNUMGEN_HPP

#include <random>
#include <chrono>

// To emulate what has been done by the professor, we make our 
// random number generator Mersenne Twister

//The usage of template, makes for a more flexible class.
//Say for example i want to generate floats, or change generator, I can easily change them now


//Bevause we're using a template class, the implementation must be in the header file
template <class ntype=double, class engine=std::mt19937_64>

class randNumGen {
    private:
        std::random_device rd;
        using myclock = std::chrono::high_resolution_clock;
        engine rng;
        static constexpr ntype m=0.0,M=1.0;
        using ud = std::uniform_real_distribution<ntype>;
        ud* unidst; // Not necessary to use a pointer here, but it is used to show how to use a pointer to uniform_real_distribution class
        public:
            void seed(int s) {
                // Usiding a seed sequence, instead of a single seed, is
                //greatly improves the quality of the random numbers generated
                //by the Mersenne Twister. seed_seq is designed to be used to take that sequence of integers and convert it into a seed sequence of "high quality"
                
                std::seed_seq seq{s+1, s+2, s+2};
                            
                //mt19937_64 method of initialization 
                rng.seed(seq);
            };
            //Clock based seed
            void rseed() {
            unsigned int t = myclock::now().time_since_epoch().count();
            std::seed_seq seq{rd() ^ t, rd() ^ t, rd() ^ t};
            rng.seed(seq);
            }

            // Generate a random number
            ntype ranf() { return (*unidst)(rng); }

            // Required RNG interface
            ntype operator()() { return (*unidst)(rng); }

            static constexpr ntype min() { return m; }
            static constexpr ntype max() { return M; }
            randNumGen() {unidst = new ud(m, M);};
            ~randNumGen() {delete unidst;};
}; 

randNumGen rng;
#endif // RANDNUMGEN_HPP