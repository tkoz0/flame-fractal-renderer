/*
Some utility functions for flame fractal renderer
*/

#pragma once

namespace tkoz
{
namespace flame
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
    return fabs(n) > bad_value_threshold<T>::value || isnan(n);
}

}
}
