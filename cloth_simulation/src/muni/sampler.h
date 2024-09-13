#pragma once
#include "common.h"

namespace muni {
struct UniformSampler {
    /** Initialize the random number generator.
            \param[in] seed The seed for the random number generator.
        */
    static void init(int seed) { srand(seed); }

    /** Generate a 1D random number in the range (0, 1).
            \return The random number.
        */
    static float next1d() { return (rand() + 0.5f) / (RAND_MAX + 1.0f); }

    /** Generate a 2D random number in the range (0, 1).
            \return The random number.
        */
    static Vec2f next2d() { return Vec2f(next1d(), next1d()); }
};
}  // namespace muni
