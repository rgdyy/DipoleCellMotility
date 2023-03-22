#ifndef DIPOLECELL_RANDOMHELPER_H
#define DIPOLECELL_RANDOMHELPER_H

#include "CommonTypes.h"
#include <random>

namespace RandomHelper {
//    enum class RandInitType : short {
//        default_init,
//        custom_seed_scalar,
//        rand_device_seed_scalar,
//        rand_device_seed_seq,
//    };

FloatingType
realUniform1D(FloatingType a,
              FloatingType b); // uniformly distributed on the interval [a, b)

IntegerType discreteUniform1D(
    IntegerType a,
    IntegerType b); // uniformly distributed on the closed interval [a, b]

bool bernoulli1D(double p);

FloatingType gaussian1D(FloatingType mean, FloatingType stddev);
} // namespace RandomHelper

#endif // DIPOLECELL_RANDOMHELPER_H
