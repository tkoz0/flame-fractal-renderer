/*
Constants for flame fractal rendering
*/

#pragma once

#include <cstdlib>

#include "types.hpp"

namespace tkoz
{
namespace flame
{

// small epsilon value for avoiding division by zero, flam3 uses 1e-10
template <typename T> struct eps {};
template <> struct eps<float>
{ static constexpr float value = 1e-10; };
template <> struct eps<double>
{ static constexpr double value = 1e-50; };

// machine epsilon (difference between 1 and the smallest greater number)
template <typename T> struct emach {};
template <> struct emach<float>  // 1.1920928955078125e-07
{ static constexpr float value = 1.0F / (float)(1  << 23); };
template <> struct emach<double> // 2.220446049250313e-16
{ static constexpr float value = 1.0 / (double)(1L << 52); };

// some limits for flame parameters
static const size_t max_dim = 65535;
template <typename T> struct max_rect {};
template <> struct max_rect<float> { static constexpr float value = 1e5F; };
template <> struct max_rect<double> { static constexpr double value = 1e10; };

// iterations to let point converge to the attractor
// paper suggests 20, using mantissa bits + 1
template <typename T> struct settle_iters {};
template <> struct settle_iters<float>
{ static constexpr size_t value = 24; };
template <> struct settle_iters<double>
{ static constexpr size_t value = 53; };

// bad value threshold, flam3 uses 1e10
template <typename T> struct bad_value_threshold {};
template <> struct bad_value_threshold<float>
{ static constexpr float value = 1e10; };
template <> struct bad_value_threshold<double>
{ static constexpr double value = 1e50; };

// multiplier to adjust numbers for scaling
template <typename T> struct scale_adjust {};
template <> struct scale_adjust<float>
{ static constexpr float value = 1.0F - 2.0F*emach<float>::value; };
template <> struct scale_adjust<double>
{ static constexpr double value = 1.0 - 2.0*emach<double>::value; };

// image scaling multiplier for [0.0,1.0] scale
// it is made to be slightly smaller than 2^bits so rounded down
// the range is [0,255] for 8 bit and [0,65535] for 16 bit
template <typename Pixel, typename Number> struct pix_scale;
template <> struct pix_scale<u8,float>
{ static constexpr float value = 256.0F * scale_adjust<float>::value; };
template <> struct pix_scale<u16,float>
{ static constexpr float value = 65536.0F * scale_adjust<float>::value; };
template <> struct pix_scale<u8,double>
{ static constexpr double value = 256.0 * scale_adjust<double>::value; };
template <> struct pix_scale<u16,double>
{ static constexpr double value = 65536.0 * scale_adjust<double>::value; };

}
}
