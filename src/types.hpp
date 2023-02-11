#pragma once

#include <cstdlib>

// forward declarations
namespace tkoz::flame
{

template <typename T, size_t N> class Affine;
template <typename num_t, size_t dims> class Flame;
template <typename num_t, size_t dims> struct IterState;
template <typename T, size_t N> class Point;
template <typename num_t, size_t dims> class XForm;

}

#include "types/affine.hpp"
#include "types/constants.hpp"
#include "types/flame.hpp"
#include "types/iter_state.hpp"
#include "types/point.hpp"
#include "types/types.hpp"
#include "types/xform.hpp"
