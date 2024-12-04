#ifndef RANDNUMGEN_HPP
#define RANDNUMGEN_HPP

#include <random>
#include <chrono>
#include <type_traits>

// Template class for random number generation
template <class ntype = double, class engine = std::mt19937_64>
class randNumGen {
private:
    using myclock = std::chrono::high_resolution_clock;
    engine rng;
    static constexpr ntype m = 0.0, M = 1.0;

    // Conditional type for distribution
    using ud = typename std::conditional<
        std::is_integral<ntype>::value,
        std::uniform_int_distribution<ntype>,
        std::uniform_real_distribution<ntype>
    >::type;

    ud unidst;

public:
    randNumGen() : unidst(m, M) {}

    void seed(int s) {
        std::seed_seq seq{s + 1, s + 2, s + 3};
        rng.seed(seq);
    }

    void rseed() {
        std::random_device rd;
        unsigned int t = myclock::now().time_since_epoch().count();
        std::seed_seq seq{rd() ^ t, rd() ^ t, rd() ^ t};
        rng.seed(seq);
    }

    ntype operator()() { return unidst(rng); }

    // Set new distribution range
    void setDistribution(ntype a, ntype b) {
        unidst = ud(a, b);
    }

    static constexpr ntype min() { return m; }
    static constexpr ntype max() { return M; }

    // Method to get a random number from the distribution
    ntype ranf() { return unidst(rng); }
};

#endif // RANDNUMGEN_HPP
