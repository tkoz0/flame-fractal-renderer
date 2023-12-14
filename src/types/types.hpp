/*
Types for convenience
*/

#pragma once

#include <cstdint>
#include <cstdlib>

#include <boost/gil.hpp>

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
using rng_t = FlameRNG<num_t,u64,4>;

// number type
// double - slower, unnecessary
// float - faster, reasonable
// (not an option but half precision may be ok for smaller buffers)
typedef double num_t;
// histogram counter type
// uint64_t - big enough for plenty
// uint32_t - reasonable choice for most renders
// uint16_t - too small to be practical (not supported)
typedef uint64_t hist_t;

static_assert(std::is_same<num_t,float>::value
           || std::is_same<num_t,double>::value);
static_assert(std::is_same<hist_t,uint32_t>::value
           || std::is_same<hist_t,uint64_t>::value);
static_assert(sizeof(num_t) == sizeof(hist_t));

// image typenames
template <typename Pixel> struct image_gray;
template <> struct image_gray<u8>
{ typedef boost::gil::gray8_image_t type; };
template <> struct image_gray<u16>
{ typedef boost::gil::gray16_image_t type; };
template <> struct image_gray<u32>
{ typedef boost::gil::gray32_image_t type; };
template <typename Pixel> struct image_rgb;
template <> struct image_rgb<u8>
{ typedef boost::gil::rgb8_image_t type; };
template <> struct image_rgb<u16>
{ typedef boost::gil::rgb16_image_t type; };
template <> struct image_rgb<u32>
{ typedef boost::gil::rgb32_image_t type; };

}
