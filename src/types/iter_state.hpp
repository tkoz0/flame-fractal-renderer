/*
Representation of the iterating point
*/

#pragma once

#include <cstdlib>

// forward declaration
namespace tkoz::flame
{
template <typename num_t, size_t dims> struct IterState;
}

#include "../settings.hpp"
#include "types.hpp"
#include "point.hpp"

namespace tkoz::flame
{

// variables must be initialized manually
template <typename num_t, size_t dims>
struct IterState
{
    static_assert(std::is_same<num_t,float>::value
               || std::is_same<num_t,double>::value);
    typedef Point<num_t,dims> point_t;
    // iterating point
    // - p = current point
    // - t = pre-affine transformation applied to p
    // - v = variation sum based on t
    point_t p,t,v;
    // rng state
    rng_t<num_t> *rng;
    // cumulative weights for xform selection
    num_t *cw;
};

#if INSTANTIATE_TEMPLATES
extern template class IterState<float,1>;
extern template class IterState<float,2>;
extern template class IterState<float,3>;
extern template class IterState<float,4>;
extern template class IterState<float,5>;

extern template class IterState<double,1>;
extern template class IterState<double,2>;
extern template class IterState<double,3>;
extern template class IterState<double,4>;
extern template class IterState<double,5>;
#endif

}
