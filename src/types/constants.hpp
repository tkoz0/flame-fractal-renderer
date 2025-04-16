/*
Constants for flame fractal rendering
*/

#pragma once

#include "types.hpp"

#include <boost/gil.hpp>

#include <cstdlib>

namespace tkoz::flame
{

// small epsilon value for avoiding division by zero, flam3 uses 1e-10
template <typename T> struct eps {};
template <> struct eps<float>
{ static constexpr float  value = 1e-10; };
template <> struct eps<double>
{ static constexpr double value = 1e-20; };
template <typename T> static constexpr auto eps_v = eps<T>::value;

// machine epsilon (difference between 1 and the smallest greater number)
template <typename T> struct emach {};
template <> struct emach<float>  // 1.1920928955078125e-07
{ static constexpr float value = 1.0F / (float)(1  << 23); };
template <> struct emach<double> // 2.220446049250313e-16
{ static constexpr float value = 1.0 / (double)(1L << 52); };
template <typename T> static constexpr auto emach_v = emach<T>::value;

// some limits for flame parameters
static const size_t max_dim = 65535;
template <typename T> struct max_rect {};
template <> struct max_rect<float>  { static constexpr float  value = 1e5F; };
template <> struct max_rect<double> { static constexpr double value = 1e10; };
template <typename T> static constexpr auto max_rect_v = max_rect<T>::value;

// iterations to let point converge to the attractor
// paper suggests 20, using mantissa bits + 1
template <typename T> struct settle_iters {};
template <> struct settle_iters<float>
{ static constexpr size_t value = 24; };
template <> struct settle_iters<double>
{ static constexpr size_t value = 53; };
template <typename T>
static constexpr auto settle_iters_v = settle_iters<T>::value;

// bad value threshold, flam3 uses 1e10
template <typename T> struct bad_value_threshold {};
template <> struct bad_value_threshold<float>
{ static constexpr float  value = 1e10; };
template <> struct bad_value_threshold<double>
{ static constexpr double value = 1e20; };
template <typename T>
static constexpr auto bad_value_threshold_v = bad_value_threshold<T>::value;

// multiplier to adjust numbers for scaling
template <typename T> struct scale_adjust_down {};
template <> struct scale_adjust_down<float>
{ static constexpr float  value = 1.0F - emach_v<float>; };
template <> struct scale_adjust_down<double>
{ static constexpr double value = 1.0 - emach_v<double>; };
template <typename T> struct scale_adjust_up {};
template <> struct scale_adjust_up<float>
{ static constexpr float  value = 1.0F + emach_v<float>; };
template <> struct scale_adjust_up<double>
{ static constexpr double value = 1.0 + emach_v<double>; };
template <typename T>
static constexpr auto scale_adjust_down_v = scale_adjust_down<T>::value;
template <typename T>
static constexpr auto scale_adjust_up_v = scale_adjust_up<T>::value;

// image scaling multiplier for [0.0,1.0] scale
// it is made to be slightly smaller than 2^bits so rounded down
// the range is [0,255] for 8 bit and [0,65535] for 16 bit
template <typename Pixel, typename Number> struct pix_scale;
template <> struct pix_scale<u8,float>
{ static constexpr float value = 256.0F * scale_adjust_down_v<float>; };
template <> struct pix_scale<u16,float>
{ static constexpr float value = 65536.0F * scale_adjust_down_v<float>; };
template <> struct pix_scale<u32,float>
{ static constexpr float value = 4294967296.0F * scale_adjust_down_v<float>; };
template <> struct pix_scale<u8,double>
{ static constexpr double value = 256.0 * scale_adjust_down_v<double>; };
template <> struct pix_scale<u16,double>
{ static constexpr double value = 65536.0 * scale_adjust_down_v<double>; };
template <> struct pix_scale<u32,double>
{ static constexpr double value = 4294967296.0 * scale_adjust_down_v<double>; };
template <typename P, typename N>
static constexpr auto pix_scale_v = pix_scale<P,N>::value;

} // namespace tkoz::flame
