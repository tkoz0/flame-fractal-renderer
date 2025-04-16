/*
Math functions
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "../types/types.hpp"

#include <ctgmath>

namespace tkoz::flame
{

// type generic sincos with function overloading

// sincos with float
static inline void sincosg(float a, float& s, float& c) { sincosf(a,&s,&c); }

// sincos with double
static inline void sincosg(double a, double& s, double& c) { sincos(a,&s,&c); }

// constants

const num_t c_pi = 3.141592653589793;
const num_t c_2pi = 6.283185307179586;
const num_t c_pi_2 = 1.5707963267948966;
const num_t c_pi_3 = 1.0471975511965976;
const num_t c_pi_4 = 0.7853981633974483;
const num_t c_pi_6 = 0.5235987755982988;

const num_t c_sqrt2 = 1.4142135623730951;
const num_t c_sqrt3 = 1.7320508075688772;
const num_t c_sqrt5 = 2.23606797749979;

const num_t c_log2 = 0.6931471805599453;
const num_t c_log3 = 1.0986122886681098;
const num_t c_log10 = 2.302585092994046;

const num_t c_e = 2.718281828459045;

} // namespace tkoz::flame
