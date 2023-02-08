/*
Representation of the iterating point
*/

#pragma once

#include <cstdlib>

#include "types.hpp"
#include "point.hpp"

namespace tkoz
{
namespace flame
{

template <typename num_t, size_t dims>
struct IterState
{
    static_assert(std::is_same<num_t,float>::value
               || std::is_same<num_t,double>::value);
    typedef Point<num_t,dims> point_t;
    // iterating point
    // - p = current point
    // - t = pre-affine transformation applied to p
    // - v = variation sum based on t
    point_t p, t, v;
    // rng state
    rng_t& rng;
    // cumulative weights for xform selection
    num_t *cw;
    // set points to 0 and use provided rng
    IterState(rng_t& rng): p(),t(),v(),rng(rng) {}
    // random true/false
    inline bool randBool()
    {
        return rng.next() & 1;
    }
    // random number in [0,1)
    inline num_t randNum()
    {
        if (sizeof(num_t) == 4) // float
            return (rng.next() >> 8) / (float)(1 << 24);
        else // double
        {
            u32 hi = rng.next() >> 6;
            u32 lo = rng.next() >> 5;
            return (((u64)hi << 27) | lo) / (double)(1LL << 53);
        }
        // code for isaac64
        //if (sizeof(num_t) == 4) // float
        //    return (rng.next() >> 40) / (float)(1 << 24);
        //else // double
        //    return (rng.next() >> 11) / (double)(1LL << 53);
    }
    // generates 2 gaussian variables (mean 0, stdev 1)
    inline void randGaussianPair(num_t& z1, num_t& z2)
    {
#if 0 // box muller transform
        num_t u1 = randNum();
        num_t u2 = (2.0*M_PI)*randNum();
        num_t r = sqrt(-2.0*log(u1));
        num_t s,c;
        sincosg(u2,&s,&c);
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
    // generates 1 gaussian variable, the unused one is wasted
    // the algorithm generates these in pairs so only use this if 1 is needed
    inline num_t randGaussian()
    {
        num_t z1,z2;
        randGaussianPair(z1,z2);
        return z1;
    }
    // random point in the biunit square/cube/hypercube
    inline point_t randPoint()
    {
        num_t x[dims];
        for (size_t i = 0; i < dims; ++i)
            x[i] = 2.0*randNum() - 1.0;
        return point_t(x);
    }
    // random point in unit square/cube/hypercube
    inline point_t randPoint2()
    {
        num_t x[dims];
        for (size_t i = 0; i < dims; ++i)
            x[i] = randNum() - 0.5;
        return point_t(x);
    }
    // random point on the unit circle/sphere/hypersphere surface
    inline point_t randDirection()
    {
        if (dims == 1)
            return point_t(copysign(1.0,randNum()-0.5));
        else if (dims == 2)
        {
            num_t a = (2.0*M_PI) * randNum();
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            return point_t(ca,sa);
        }
        else if (dims == 3)
        {
            num_t p = acos(2.0*randNum()-1.0);
            num_t t = (2.0*M_PI) * randNum();
            num_t st,ct,sp,cp;
            sincosg(p,&sp,&cp);
            sincosg(t,&st,&ct);
            return point_t(st*cp,st*sp,ct);
        }
        else
        {
            num_t x[dims];
            for (size_t i = 0; i+1 < dims; i += 2)
                randGaussianPair(x[i],x[i+1]);
            if (dims % 2)
                x[dims-1] = randGaussian();
            point_t p(x);
            return p/p.norm2();
        }
    }
    // use cumulative weights to select xform
    inline u32 randXFormIndex()
    {
        u32 ret = 0;
        num_t r = randNum();
        while (cw[ret] < r)
            ++ret;
        return ret;
    }
};

}
}
