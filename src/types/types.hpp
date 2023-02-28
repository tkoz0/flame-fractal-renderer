/*
Types for convenience
*/

#pragma once

#include <cstdint>
#include <cstdlib>

#include "../rng/flame_rng.hpp"

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

// random number generator type to use
// rparam=4 recommended for simulations, rparam=8 recommended for cryptography
// state size is 2^rparam of word_t
// word_t can be u32 for isaac32 or u64 for isaac64
template <typename num_t>
using rng_t = FlameRNG<num_t,u32,4>;

}
