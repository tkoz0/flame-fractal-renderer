/*
Types for convenience
*/

#pragma once

#include <cstdint>
#include <cstdlib>

namespace tkoz::flame
{

// integer types
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;

// number type
// double - slightly slower, unnecessary
// float - slightly faster, reasonable
// (not an option but half precision may be ok for smaller buffers)
typedef double num_t;
// histogram counter type
// uint64_t - big enough for plenty
// uint32_t - reasonable choice for most renders
// uint16_t - too small to be practical (not supported)
typedef uint64_t hist_t;

}

#include "../rng/flame_rng.hpp"

namespace tkoz::flame
{

static_assert(std::is_same<num_t,float>::value
           || std::is_same<num_t,double>::value);
static_assert(std::is_same<hist_t,uint32_t>::value
           || std::is_same<hist_t,uint64_t>::value);
static_assert(sizeof(num_t) == sizeof(hist_t));

// random number generator type to use
// rparam=4 recommended for simulations, rparam=8 recommended for cryptography
// state size is 2^rparam of word_t
// word_t can be u32 for isaac32 or u64 for isaac64
using rng_t = FlameRNG<num_t,hist_t,4>;

}
