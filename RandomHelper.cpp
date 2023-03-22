#include "RandomHelper.h"
#include "pcg-cpp-0.98/pcg_random.hpp"
#include <array>
#include <fstream>
#include <functional>
#include <iostream>

namespace RandomHelper {
//    auto newRandSeedOfstream() {
//        std::ofstream rand_ostream("rand_seed.txt");
//        rand_ostream.exceptions(~std::ofstream::goodbit);
//        return rand_ostream;
//    }

auto newURNG_defaultInit() {
    std::default_random_engine random_engine;
    std::cout << "Random Engine default initialized.\n";
    return random_engine;
}

auto newURNG_customSeedInit(unsigned long seed) {
    std::default_random_engine random_engine(seed);
    //        std::cout << "Random Engine seed = " << seed << '\n';
    return random_engine;
}

auto newURNG_randomDeviceInit() {
    auto seed = std::random_device()();
    std::default_random_engine random_engine(seed);
    //        std::cout << "Random Engine seed = " << seed << '\n';
    return random_engine;
}

auto newURNG_mt19937SeedSeqInit() {
    std::random_device random_device;
    std::array<std::random_device::result_type, std::mt19937::state_size>
        seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(random_device));
    std::seed_seq seed_seq(seed_data.begin(), seed_data.end());
    std::mt19937 random_engine(seed_seq);
    //        auto &&rand_seed_ostream = std::cout;
    //        for (auto s : seed_data) {
    //            rand_seed_ostream << s << '\n';
    //        }
    return random_engine;
}

auto newURNG_pcgDefaultInit() {
    //        https://www.pcg-random.org/using-pcg-cpp.html
    return pcg32(); // Make a random number engine
}

auto newURNG_pcgSeedSeqInit() {
    //        https://www.pcg-random.org/using-pcg-cpp.html
    pcg_extras::seed_seq_from<std::random_device>
        seed_source;           // Seed with a real random value, if available
    return pcg32(seed_source); // Make a random number engine
}

auto global_random_engine = newURNG_pcgDefaultInit();

FloatingType realUniform1D(FloatingType a, FloatingType b) {
    static std::uniform_real_distribution<FloatingType>
        uniform_real_distribution;
    return uniform_real_distribution(
        global_random_engine,
        decltype(uniform_real_distribution)::param_type(a, b));
}

IntegerType discreteUniform1D(IntegerType a, IntegerType b) {
    static std::uniform_int_distribution<IntegerType> uniform_int_distribution;
    return uniform_int_distribution(
        global_random_engine,
        decltype(uniform_int_distribution)::param_type(a, b));
}

bool bernoulli1D(double p) {
    static std::bernoulli_distribution bernoulli_distribution;
    return bernoulli_distribution(
        global_random_engine, decltype(bernoulli_distribution)::param_type(p));
}

FloatingType gaussian1D(FloatingType mean, FloatingType stddev) {
    static std::normal_distribution<FloatingType> normal_distribution;
    return normal_distribution(
        global_random_engine,
        decltype(normal_distribution)::param_type(mean, stddev));
}
} // namespace RandomHelper
