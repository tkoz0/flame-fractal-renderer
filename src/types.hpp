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

// macros for using SFINAE
#define ENABLE_IF(COND,RET) template <typename RET2 = RET> \
    typename std::enable_if<(COND),RET2>::type

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

// variadic template all same type
// https://www.fluentcpp.com/2019/01/25/variadic-number-function-parameters-type/
template <bool...> struct bool_pack{};
template <typename...Ts>
using vconj = std::is_same<bool_pack<true,Ts::value...>,
                           bool_pack<Ts::value...,true>>;
// all of Ts same as T
template <typename T, typename...Ts>
using vallsame = vconj<std::is_same<T,Ts>...>;
// all of Ts convertible to T
template <typename T, typename...Ts>
using vallconv = vconj<std::is_convertible<T,Ts>...>;

// similar to std::enable_if but for > 2 cases of function specialization
template <size_t N1, size_t N2, typename T = void>
struct enable_if_eq {};
template <size_t N, typename T>
struct enable_if_eq<N,N,T> { typedef T type; };
#define ENABLE_IFEQ(N1,N2,RET) template <typename RET2 = RET> \
    typename enable_if_eq<N1,N2,RET2>::type

// point in N dimension space, using number type T
template <typename T, size_t N>
class Point
{
    static_assert(N > 0 && N < 100);
private:
    std::array<T,N> vec;
public:
    inline Point()
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = 0;
    }
    inline Point(const T x_[N]) // from length N array
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = x_[i];
    }
    inline Point(const Json& j)
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
    inline Point(const std::array<T,N>& x_): vec(x_) {}
    // Variadic templated constructor with N arguments
    // U,Us... to require >= 1 argument so no ambiguity results
    // argumuent types can be anything that can be casted to type T
    // using dummy argumuent for SFINAE
    // - U (1) + Us (N-1) == length N
    // - U and Us... are convertible to T
    template <typename U, typename...Us, typename std::enable_if<
        sizeof...(Us) == N-1
        && std::is_convertible<U,T>::value
        && vallconv<T,Us...>::value, size_t>::type = 0>
    inline Point(U u, Us... us)
    {
        // explicitly cast all to type T for initializer list
        vec = {(T) u, (T) us...};
    }
    inline Point(const Point<T,N>& p)
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = p.vec[i];
    }
    inline Point<T,N>& operator=(const Point<T,N>& p)
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = p.vec[i];
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
    // x value (1st component)
    ENABLE_IF(N>=1,T) inline x() const
    {
        return vec[0];
    }
    // y value (2nd component) (requires 2D)
    ENABLE_IF(N>=2,T) inline y() const
    {
        return vec[1];
    }
    // z value (3rd component) (requires 3D)
    ENABLE_IF(N>=3,T) inline z() const
    {
        return vec[2];
    }
    ENABLE_IF(N>=2,void) inline getXY(T& x_, T& y_) const
    {
        x_ = x();
        y_ = y();
    }
    ENABLE_IF(N>=3,void) inline getXYZ(T& x_, T& y_, T& z_) const
    {
        getXY(x_,y_);
        z_ = z();
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
    // fold left
    inline T foldl(std::function<T(T,T)> func, T init = 0) const
    {
        for (size_t i = 0; i < N; ++i)
            init = func(init,vec[i]);
        return init;
    }
    // fold right
    inline T foldr(std::function<T(T,T)> func, T init = 0) const
    {
        for (size_t i = N-1; i--;)
            init = func(vec[i],init);
        return init;
    }
    // 2-norm functions
    ENABLE_IFEQ(N,1,T) inline norm2() const { return abs(x()); }
    ENABLE_IFEQ(N,2,T) inline norm2() const { return hypot(x(),y()); }
    ENABLE_IFEQ(N,3,T) inline norm2() const { return sqrt(norm2sq()); }
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
            ret += abs(vec[i]);
        return ret;
    }
    // inf-norm, TODO ensure optimized for small N
    inline T norminf() const
    {
        T ret = 0;
        for (size_t i = 0; i < N; ++i)
            ret = std::max(ret,abs(vec[i]));
        return ret;
    }
    // angle in xy plane
    ENABLE_IF(N==2||N==3,T) inline angle() const
    {
        return atan2(y(),x());
    }
    // angle of inclination (with +z axis) in 3d sphereical coordinates
    ENABLE_IF(N==3,T) inline inclination() const
    {
        return acos(z()/norm2());
    }
    // polar coordinates
    ENABLE_IF(N==2,void) inline polar(T& radius, T& theta) const
    {
        radius = norm2();
        theta = angle();
    }
    // spherical coordinates
    ENABLE_IF(N==3,void) inline spherical(T& radius, T& theta, T& phi) const
    {
        radius = norm2();
        theta = inclination();
        phi = angle();
    }
    // get radius and sin,cos of angle
    ENABLE_IF(N==2,void) inline getRadiusSinCos(T& r, T& s, T& c) const
    {
        r = norm2();
        s = y() / r;
        c = x() / r;
    }
    // get sin,cos of angle
    ENABLE_IF(N==2,void) inline getSinCos(T& s, T& c) const
    {
        T r;
        getRadiusSinCos(r,s,c);
    }
};

// affine transformation in N dimensions
template <typename T, size_t N>
class Affine
{
    static_assert(N > 0 && N < 100);
private:
    std::array<Point<T,N>,N> A; // linear transformation A*x + b
    Point<T,N> b;               // A as array of its rows
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
    inline const std::array<Point<T,N>,N>& getA() const
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

#undef ENABLE_IF
