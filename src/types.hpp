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

template <typename T> inline bool bad_value(T n)
{
    return fabs(n) > bad_value_threshold<T>::value || isnan(n);
}

template <typename T, size_t N> struct Affine;
template <typename T, size_t N> struct Point;

// point in N dimension space, using number type T
template <typename T, size_t N>
class Point
{
    static_assert(N > 0 && N < 100);
private:
    T vec[N];
public:
    Point()
    {
        memset(vec,0,sizeof(vec));
    }
    Point(const T *x_) // from length N array
    {
        memcpy(vec,x_,sizeof(vec));
    }
    // TODO remove, temporary to make variations.hpp compile
    /*
    Point(T x, T y) { vec[0] = x; vec[1] = y; }
    inline T r() { return hypot(vec[0],vec[1]); }
    inline T r2() { return vec[0]*vec[0] + vec[1]*vec[1]; }
    inline T atanxy() { return atan2(vec[0],vec[1]); }
    inline T atanyx() { return atan2(vec[1],vec[0]); }
    */
    Point(const Json& j)
    {
        JsonArray a = j.arrayValue();
        if (a.size() != N)
            throw std::runtime_error("point: incorrect size");
        for (size_t i = 0; i < N; ++i)
        {
            if (a[i].isInt())
                vec[i] = (T) a[i].intValue();
            else if (a[i].isFloat())
                vec[i] = (T) a[i].floatValue();
            else
                throw std::runtime_error("point: not a number");
        }
    }
    // TODO constructor taking N args of type T
    Point(const Point<T,N>& p)
    {
        memcpy(vec,p.vec,sizeof(vec));
    }
    inline Point<T,N>& operator=(const Point<T,N>& p)
    {
        memcpy(vec,p.vec,sizeof(vec));
        return *this;
    }
    inline Point<T,N>& operator+=(const Point<T,N>& p)
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] += p.vec[i];
        return *this;
    }
    inline Point<T,N>& operator-=(const Point<T,N>& p)
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] -= p.vec[i];
        return *this;
    }
    inline Point<T,N>& operator*=(T k)
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] *= k;
        return *this;
    }
    inline Point<T,N>& operator/=(T k)
    {
        return *this *= (1/k);
    }
    inline T operator[](size_t i) const
    {
        return vec[i];
    }
    inline T& operator[](size_t i)
    {
        return vec[i];
    }
    inline T x() const
    {
        return vec[0];
    }
    inline T y() const
    {
        return vec[1];
    }
    inline T z() const
    {
        return vec[2];
    }
    friend inline Point<T,N> operator+(const Point<T,N>& a,
                                       const Point<T,N>& b)
    {
        Point<T,N> ret = a;
        ret -= b;
        return ret;
    }
    friend inline Point<T,N> operator-(const Point<T,N>& a,
                                       const Point<T,N>& b)
    {
        Point<T,N> ret = a;
        ret -= b;
        return ret;
    }
    friend inline Point<T,N> operator*(const Point<T,N>& a, T k)
    {
        Point<T,N> ret = a;
        ret *= k;
        return ret;
    }
    friend inline Point<T,N> operator*(T k, const Point<T,N>& a)
    {
        return a*k;
    }
    // dot product
    friend inline T operator*(const Point<T,N>& a, const Point<T,N>& b)
    {
        T ret = 0;
        for (size_t i = 0; i < N; ++i)
            ret += a[i] * b[i];
        return ret;
    }
    friend inline Point<T,N> operator/(const Point<T,N>& a, T k)
    {
        Point<T,N> ret = a;
        ret /= k;
        return ret;
    }
    friend inline Point<T,N> operator/(T k, const Point<T,N>& a)
    {
        return a/k;
    }
    // map each coordinate under the same function
    inline Point<T,N> map(std::function<T(T)> func) const
    {
        Point<T,N> ret;
        for (size_t i = 0; i < N; ++i)
            ret[i] = func(vec[i]);
        return ret;
    }
    // 2-norm, TODO hypot for N=2
    inline T norm2() const
    {
        return sqrt(norm2sq());
    }
    // 2-norm squared, TODO ensure loop is unrolled for small N
    inline T norm2sq() const
    {
        T ret = 0;
        for (size_t i = 0; i < N; ++i)
            ret += vec[i]*vec[i];
        return ret;
    }
    // 1-norm, TODO ensure optimized for small N
    inline T norm1() const
    {
        T ret = 0;
        for (size_t i = 0; i < N; ++i)
            ret += vec[i];
        return ret;
    }
    // inf-norm, TODO ensure optimized for small N
    inline T norminf() const
    {
        T ret = 0;
        for (size_t i = 0; i < N; ++i)
            ret = std::max(ret,vec[i]);
        return ret;
    }
};

// affine transformation in N dimensions
template <typename T, size_t N>
class Affine
{
    static_assert(N > 0 && N < 100);
private:
    Point<T,N> A[N]; // linear transformation A*x + b
    Point<T,N> b;    // A as array of its rows
public:
    Affine(): A(), b()
    {
        for (size_t i = 0; i < N; ++i)
            A[i][i] = 1;
    }
    Affine(const T *A_, const T *b_) // from length N*N and N arrays
    {
        b = Point<T,N>(b_);
        for (size_t i = 0; i < N; ++i)
            A[i] = Point<T,N>(A_+(N*i));
    }
    Affine(const Json& j) // {"A": [2d array], "b": [1d array]}}
    {
        JsonObject jo = j.objectValue();
        if (jo.find("A") == jo.end()) // use identity matrix
        {
            for (size_t i = 0; i < N; ++i)
            {
                A[i] = Point<T,N>();
                A[i][i] = 1;
            }
        }
        else
        {
            JsonArray ja = jo["A"].arrayValue();
            if (ja.size() != N)
                throw std::runtime_error("affine: incorrect size");
            for (size_t i = 0; i < N; ++i)
                A[i] = Point<T,N>(ja[i]);
        }
        if (jo.find("b") == jo.end()) // use 0
            b = Point<T,N>();
        else
            b = Point<T,N>(jo["b"]);
    }
    inline const Point<T,N>& operator[](size_t i) const
    {
        return A[i];
    }
    inline const Point<T,N>* getA() const
    {
        return A;
    }
    inline const Point<T,N>& getB() const
    {
        return b;
    }
    // TODO ensure loop is unrolled for small N
    inline Point<T,N> apply_to(const Point<T,N>& x) const
    {
        Point<T,N> ret(b);
        for (size_t i = 0; i < N; ++i)
            ret[i] += A[i] * x;
        return ret;
    }
};

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Point<T,N>& p)
{
    os << "(" << p[0];
    for (size_t i = 1; i < N; ++i)
        os << "," << p[i];
    os << ")";
    return os;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Affine<T,N>& a)
{
    os << "{A=[" << a[0];
    for (size_t i = 1; i < N; ++i)
        os << "," << a[i];
    os << "],b=" << a.getB() << "}";
    return os;
}

// forward declarations for renderer.hpp
template <typename num_t, size_t dims, typename rand_t> struct XFormVar;
template <typename num_t, size_t dims, typename rand_t> class XForm;
template <typename num_t, size_t dims, typename rand_t> class Flame;
template <typename num_t, typename hist_t, typename rand_t>
class RendererBasic;

// iterator state used by variation functions, defined with renderer
template <typename num_t,
          size_t dims = 2,
          typename rand_t = Isaac<u32,4>> struct IterState;

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
