#pragma once

#include <cstdlib>

// forward declarations
namespace tkoz::flame::vars
{

template <typename num_t, size_t dims> class Variation;

}

// base variation class
#include "variations/variation_base.hpp"

// individual variations
#include "variations/variations.hpp"

// factory depends on the individual variations
// so define its implementation at the end
#include "variations/variation_factory.hpp"
