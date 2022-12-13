/*
Types for flame fractal renderer. In templates, num_t is float or double.
*/

#pragma once

#include <cstdint>
#include <ctgmath>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "isaac.hpp"
#include "jrand.hpp"
#include "json_small.hpp"

namespace tkoz
{
namespace flame
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

// type generic sincos for float and double
static inline void sincosg(float a, float *s, float *c) { sincosf(a,s,c); }
static inline void sincosg(double a, double *s, double *c) { sincos(a,s,c); }

// some limits for flame parameters
static const size_t max_dim = 64000;
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
{ static constexpr float value = 1.0F - 8.0F*emach<float>::value; };
template <> struct scale_adjust<double>
{ static constexpr double value = 1.0 - 8.0*emach<double>::value; };

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

template <typename T> inline bool bad_value(T n)
{
    return fabs(n) > bad_value_threshold<T>::value || isnan(n);
}

template <typename T> struct rgb_t
{
    T r, g, b;
    rgb_t(): r(0),g(0),b(0) {}
    rgb_t(T r, T g, T b): r(r),g(g),b(b) {}
};

template <typename T> struct rgba_t
{
    T r, g, b, a;
    rgba_t(): r(0),g(0),b(b),a(0) {}
    rgba_t(T r, T g, T b, T a): r(r),g(g),b(b),a(a) {}
};

template <typename T> struct Affine2D;

// point in 2d space
template <typename T> struct Point2D
{
    T x, y;
    Point2D(): x(0),y(0) {}
    Point2D(T x, T y): x(x),y(y) {}
    inline Point2D& operator+=(const Point2D& p)
    { x += p.x; y += p.y; return *this; }
    inline Point2D& operator-=(const Point2D& p)
    { x -= p.x; y -= p.y; return *this; }
    inline Point2D& operator*=(T k)
    { x *= k; y *= k; return *this; }
    inline Point2D& operator/=(T k)
    { x /= k; y /= k; return *this; }
    inline T r() const { return hypot(x,y); }
    inline T r2() const { return x*x + y*y; }
    inline T atan() const { return atan2(x,y); }
    inline T atanyx() const { return atan2(y,x); }
    // modify this point by an affine transformation
    inline Point2D<T> transform(const Affine2D<T>& t)
    {
        T xnew = t.a*x+t.b*y+t.c;
        T ynew = t.d*x+t.e*y+t.f;
        x = xnew;
        y = ynew;
    }
    friend inline Point2D<T> operator+(Point2D<T> a, Point2D<T> b)
    { Point2D<T> ret = a; ret += b; return ret; }
    friend inline Point2D<T> operator-(Point2D<T> a, Point2D<T> b)
    { Point2D<T> ret = a; ret -= b; return ret; }
    friend inline Point2D<T> operator*(Point2D<T> a, T k)
    { Point2D<T> ret = a; ret *= k; return ret; }
    friend inline Point2D<T> operator*(T k, Point2D<T> a)
    { return a*k; }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Point2D<T>& p)
{
    os << "(" << p.x << "," << p.y << ")";
    return os;
}

// parameters for affine transformation (x,y) -> (a*x+b*y+c,d*x+e*y+f)
template <typename T> struct Affine2D
{
    T a, b, c, d, e, f;
    // default to (x,y) -> (x,y)
    Affine2D(): a(1.0),b(0.0),c(0.0),d(0.0),e(1.0),f(0.0) {}
    // take 6 arguments
    Affine2D(T a, T b, T c, T d, T e, T f):
        a(a),b(b),c(c),d(d),e(e),f(f) {}
    inline Point2D<T> apply_to(Point2D<T> p) const
    { return Point2D<T>(a*p.x+b*p.y+c, d*p.x+e*p.y+f); }
    inline Point2D<T> apply_to(T x, T y) const
    { return Point2D<T>(a*x+b*y+c, d*x+e*y+f); }
};

// forward declarations for renderer.hpp
template <typename num_t, typename rand_t> struct XFormVar;
template <typename num_t, typename rand_t> class XForm;
template <typename num_t, typename rand_t> class Flame;
template <typename num_t, typename hist_t, typename rand_t>
class RendererBasic;

// precalc flags for renderer and variations
/*
unused code from previous C version
#define PC_THETA (1 << 0)
#define PC_PHI   (1 << 1)
#define PC_SINT  (1 << 2)
#define PC_COST  (1 << 3)
#define PC_R     (1 << 4)
#define PC_R2    (1 << 5)
*/

}
}
