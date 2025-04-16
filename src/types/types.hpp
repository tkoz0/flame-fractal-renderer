/*
Types for convenience
*/

#pragma once

#include <cstdint>
#include <cstdlib>
#include <type_traits>

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

static_assert(std::is_same_v<num_t,float>
           || std::is_same_v<num_t,double>);
static_assert(std::is_same_v<hist_t,uint32_t>
           || std::is_same_v<hist_t,uint64_t>);
static_assert(sizeof(num_t) == sizeof(hist_t));

} // namespace tkoz::flame
