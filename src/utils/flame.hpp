/*
Some utility functions for flame fractal renderer
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

#include "../types/constants.hpp"

namespace tkoz::flame
{

// max integer value stored in another integer type
template <typename T, typename U> struct max_int_as
{
    static_assert(sizeof(T) < sizeof(U));
    static constexpr U value = 1uLL << (8*sizeof(T));
};

// bad value is NaN or absolute value above the threshold
template <typename T> inline bool bad_value(T n)
{
    return std::fabs(n) > bad_value_threshold_v<T> || std::isnan(n);
}

} // namespace tkoz::flame
