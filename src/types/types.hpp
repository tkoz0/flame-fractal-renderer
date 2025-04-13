/*
Types for convenience
*/

#pragma once

#include <cstdint>
#include <cstdlib>

namespace tkoz::flame
{

// integer types
using u8  = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

// number type
// double - slightly slower, unnecessary
// float - slightly faster, reasonable
// (not an option but half precision may be ok for smaller buffers)
using num_t = double;

// histogram counter type
// uint64_t - big enough for plenty
// uint32_t - reasonable choice for most renders
// uint16_t - too small to be practical (not supported)
using hist_t = uint64_t;

}

#include "../rng/flame_rng.hpp"

namespace tkoz::flame
{

static_assert(std::is_same_v<num_t,float>
           || std::is_same_v<num_t,double>);
static_assert(std::is_same_v<hist_t,uint32_t>
           || std::is_same_v<hist_t,uint64_t>);
static_assert(sizeof(num_t) == sizeof(hist_t));

// random number generator type to use
// rparam=4 recommended for simulations, rparam=8 recommended for cryptography
// state size is 2^rparam of word_t
// word_t can be u32 for isaac32 or u64 for isaac64
using rng_t = FlameRNG<num_t,hist_t,4>;

} // namespace tkoz::flame
