/*
Representation of a point in N dimensional real space.
*/

#pragma once

#include <cstdlib>
#include <iostream>

#include "../utils/sfinae.hpp"

// macros for using SFINAE
#define ENABLE_IF(COND,RET) template <typename RET2 = RET> \
    typename std::enable_if<(COND),RET2>::type

namespace tkoz
{
namespace flame
{

// point in N dimension space, using number type T
template <typename T, size_t N>
class Point
{
    static_assert(N > 0 && N < 100);
private:
    std::array<T,N> vec;
    typedef Point<T,N> self_type; // for ENABLE_IF macro
public:
    /* === constructors === */
    // initialize the point to the zero vector
    inline Point()
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = 0;
    }
    // initialize from length N array
    inline Point(const T x_[N])
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = x_[i];
    }
    // initialize from JSON array
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
    // initialize from std::array object
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
    // copy constructor
    inline Point(const Point<T,N>& p)
    {
        for (size_t i = 0; i < N; ++i)
            vec[i] = p.vec[i];
    }
    /* === assignment operators === */
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
    /* === element access === */
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
    // set first 2 components with references
    ENABLE_IF(N>=2,void) inline getXY(T& x_, T& y_) const
    {
        x_ = x();
        y_ = y();
    }
    // set first 3 components with references
    ENABLE_IF(N>=3,void) inline getXYZ(T& x_, T& y_, T& z_) const
    {
        getXY(x_,y_);
        z_ = z();
    }
    // return the ith component (0 indexed)
    inline T get(size_t i) const
    {
        return vec[i];
    }
    // return interval array
    inline const std::array<T,N>& getArray() const
    {
        return vec;
    }
    /* === arithmetic operators === */
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
    /* === functional programming stuff === */
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
    /* === vector norm functions === */
    // vector p-norms (p == 0 means infinity norm)
    template <u32 p> inline T norm() const
    {
        T ret;
        if (N == 1)
            return abs(x());
        switch (p)
        {
        case 0: // infinity norm
            ret = abs(x());
            for (size_t i = 1; i < N; ++i)
                ret = std::max(ret,abs(vec[i]));
            return ret;
        case 1: // l1 norm
            return normsum<p>();
        case 2: // l2 norm
            return sqrt(normsum<p>());
        default: // p-norm
            return pow(normsum<p>(),1/p);
        }
    }
    inline T norm(T p) const
    {
        return pow(normsum(p),1/p);
    }
    // vector p-norms before taking the root (p != 0)
    template <u32 p> inline T normsum() const
    {
        T ret;
        switch (p)
        {
        case 1: // l1 norm
            ret = abs(x());
            for (size_t i = 1; i < N; ++i)
                ret += abs(vec[i]);
            return ret;
        case 2: // l2 norm
            ret = x()*x();
            for (size_t i = 1; i < N; ++i)
                ret += vec[i]*vec[i];
            return ret;
        default: // p-norm
            ret = pow(abs(x()),(T)p);
            for (size_t i = 1; i < N; ++i)
                ret += pow(abs(vec[i]),(T)p);
            return ret;
        }
    }
    inline T normsum(T p) const
    {
        T ret = pow(abs(x()),p);
        for (size_t i = 1; i < N; ++i)
            ret += pow(abs(vec[i]),p);
        return ret;
    }
    // 1-norm
    inline T norm1() const { return norm<1>(); }
    // 2-norm
    inline T norm2() const { return norm<2>(); }
    // infinity-norm
    inline T norminf() const { return norm<0>(); }
    inline T normmax() const { return norm<0>(); }
    // 2-norm squared
    inline T norm2sq() const { return normsum<2>(); }
    /* === coordinate systems === */
    // angle in a plane [0,2pi)
    inline T angle(size_t d1 = 0, size_t d2 = 1) const
    {
        return atan2(vec[d2],vec[d1]);
    }
    // alternative angle in a plane [-pi,pi]
    inline T angle2(size_t d1 = 0, size_t d2 = 1) const
    {
        bool sign = std::signbit(vec[d2]);
        static const T signtable[2] = {1,-1};
        return signtable[sign] * acos(vec[d1]/hypot(vec[d1],vec[d2]));
    }
    // polar coordinates in 2d (standard polar system)
    ENABLE_IF(N==2,void) inline getPolar(T& radius, T& theta) const
    {
        radius = norm2();
        theta = angle(); // [0,2pi)
    }
    // spherical coordinates in 3d (system commonly used in physics)
    ENABLE_IF(N==3,void) inline getSpherical(T& radius, T& theta, T& phi) const
    {
        radius = norm2();
        theta = acos(z()/radius); // [0,pi]
        phi = angle(); // [0,2pi)
    }
    // generalized spherical coordinates (see wikipedia, x_i order reversed)
    // ret[0] = radius, ret[1..N-2] in [0,pi], ret[N-1] in [-pi,pi]
    ENABLE_IF(N>=2,self_type) inline toSpherical() const
    {
        Point<T,N> ret;
        T sqsum = vec[0]*vec[0] + vec[1]*vec[1];
        static const T signtable[2] = {1,-1};
        bool sign = std::signbit(vec[1]);
        ret[N-1] = signtable[sign] * acos(vec[1]/sqrt(sqsum)); // [-pi,pi]
        for (size_t i = 2; i < N; ++i)
        {
            sqsum += vec[i]*vec[i];
            ret[N-i] = acos(vec[i]/sqrt(sqsum)); // [0,pi]
        }
        ret[0] = sqrt(sqsum);
        return ret;
    }
    // set this point based on spherical coordinates
    ENABLE_IF(N>=2,void) inline setFromSpherical(const Point<T,N>& c)
    {
        T rsinprod = c[0]; // r
        vec[N-1] = rsinprod * cos(c[1]);
        for (size_t i = 1; i < N-1; ++i)
        {
            rsinprod *= sin(c[i]);
            vec[N-1-i] = rsinprod * cos(c[i+1]);
        }
        vec[0] = rsinprod * sin(c[N-1]);
    }
    /* === calculate parameters for use in variation functions === */
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

// multiply components similar to dot product but not adding them
template <typename T, size_t N>
Point<T,N> multComponents(const Point<T,N> a, const Point<T,N> b)
{
    Point<T,N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = a[i] * b[i];
    return ret;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Point<T,N>& p)
{
    os << "(" << p[0];
    for (size_t i = 1; i < N; ++i)
        os << "," << p[i];
    os << ")";
    return os;
}

}
}
