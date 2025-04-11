/*
Random number generator to use for flame fractal iteration
Based on ISAAC with extra helpful functions
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

#include "../utils/math.hpp"
#include "../utils/sfinae.hpp"

// forward declaration
namespace tkoz::flame
{
template <typename num_t, typename word_t, size_t rparam>
class FlameRNG;
}

#include "isaac.hpp"
#include "../types/point.hpp"
#include "../types/types.hpp"
#include "../utils/sfinae.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz::flame
{

template <typename num_t, typename word_t, size_t rparam>
class FlameRNG: public Isaac<word_t,rparam>
{
    static_assert(std::is_same_v<num_t,float>
               || std::is_same_v<num_t,double>);

public:

    FlameRNG(): Isaac<word_t,rparam>() {}

    FlameRNG(u32 seed): Isaac<word_t,rparam>(seed) {}

    FlameRNG(u64 seed): Isaac<word_t,rparam>(seed) {}

    inline word_t nextWord()
    {
        return Isaac<word_t,rparam>::next();
    }

    // random true/false
    inline bool randBool()
    {
        return nextWord() & 1;
    }

    // random number in [0,1)
    inline num_t randNum()
    {
        if (std::is_same_v<num_t,float>)
        {
            if (std::is_same_v<word_t,u32>)
                return (nextWord() >> 8) / (float)(1 << 24);
            else // u64
                return (nextWord() >> 40) / (float)(1 << 24);
        }
        else // double
        {
            if (std::is_same_v<word_t,u32>)
            {
                u32 hi = nextWord() >> 6;
                u32 lo = nextWord() >> 5;
                return (((u64)hi << 27) | lo) / (double)(1LL << 53);
            }
            else // u64
                return (nextWord() >> 11) / (double)(1LL << 53);
        }
    }

    // random integer in [0,n)
    // uniformly distributed only when n is a power of 2
    inline word_t randInt(word_t n)
    {
        return nextWord() % n;
    }

    // generates 2 gaussian variables
    inline void randGaussianPair(num_t& z1, num_t& z2)
    {
#if 1 // box muller transform
        num_t u1 = randNum();
        // TODO generate u2 by dividing the integer by approx max_int/2pi
        // or is it faster to divide by the power of 2 and then multiply 2pi
        num_t u2 = (2.0*M_PI)*randNum();
        num_t r = sqrt(-2.0*log(u1));
        num_t s,c;
        sincosg(u2,s,c);
        z1 = r*c;
        z2 = r*s;
#else // marsaglia polar method
        num_t v1,v2,s,r;
        do
        {
            v1 = 2.0*randNum() - 1.0;
            v2 = 2.0*randNum() - 1.0;
            s = v1*v1 + v2*v2;
        }
        while (s >= 1.0);
        r = sqrt(-2.0*log(s)/s);
        z1 = v1*r;
        z2 = v2*r;
#endif
    }

    // generates 1 gaussian variable, the unused 1 is wasted
    inline num_t randGaussian()
    {
        num_t z1,z2;
        randGaussianPair(z1,z2);
        return z1;
    }

    // random point in [-1,1]x[-1,1]x...
    template <size_t dims>
    inline Point<num_t,dims> randPoint()
    {
        Point<num_t,dims> ret;
        for (size_t i = 0; i < dims; ++i)
            ret[i] = 2.0*randNum() - 1.0;
        return ret;
    }

    // random point in [-0.5,0.5]x[-0.5,0.5]x...
    template <size_t dims>
    inline Point<num_t,dims> randPoint2()
    {
        Point<num_t,dims> ret;
        for (size_t i = 0; i < dims; ++i)
            ret[i] = randNum() - 0.5;
        return ret;
    }

    // random point on the unit 0-sphere (random sign)
    template <size_t dims> inline
    enable_if_eq_t<dims,1,Point<num_t,1>> randDirection()
    {
        return Point<num_t,1>(copysign(1.0,randNum()-0.5));
    }

    // random point on the unit 1-sphere (circle)
    template <size_t dims> inline
    enable_if_eq_t<dims,2,Point<num_t,2>> randDirection()
    {
        num_t a = (2.0*M_PI) * randNum();
        num_t sa,ca;
        sincosg(a,sa,ca);
        return Point<num_t,2>(ca,sa);
    }

    // random point on the unit 2-sphere (sphere)
    template <size_t dims> inline
    enable_if_eq_t<dims,3,Point<num_t,3>> randDirection()
    {
#if 0 // first method on wolfram mathworld
        num_t p = acos(2.0*randNum() - 1.0);
        num_t t = (2.0*M_PI) * randNum();
        num_t st,ct,sp,cp;
        sincosg(p,sp,cp);
        sincosg(t,st,ct);
        return Point<num_t,3>(st*cp,st*sp,ct);
#else // second method on wolfram mathworld
        num_t u = 2.0*randNum() - 1.0;
        num_t t = (2.0*M_PI) * randNum();
        num_t r = sqrt(1.0 - u*u);
        num_t st,ct;
        sincosg(t,st,ct);
        return Point<num_t,3>(r*ct,r*st,u);
#endif
    }

    // random point on the unit (dims-1)-sphere
    template <size_t dims> inline
    std::enable_if_t<(dims>3),Point<num_t,dims>> randDirection()
    {
        // gaussian distributed points and then normalize
        Point<num_t,dims> ret;
        for (size_t i = 0; i+1 < dims; i += 2)
            randGaussianPair(ret[i],ret[i+1]);
        if (dims % 2)
            ret[dims-1] = randGaussian();
        return ret/ret.norm2();
    }
};

} // namespace tkoz::flame

#undef likely
#undef unlikely
