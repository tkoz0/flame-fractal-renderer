/*
Math functions
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

namespace tkoz
{
namespace flame
{

// type generic sincos with function overloading

// sincos with float
static inline void sincosg(float a, float *s, float *c) { sincosf(a,s,c); }

// sincos with double
static inline void sincosg(double a, double *s, double *c) { sincos(a,s,c); }

}
}
