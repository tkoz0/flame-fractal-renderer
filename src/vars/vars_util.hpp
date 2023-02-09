/*
Helper functions for dealing with variations
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

#include "../types/types.hpp"
#include "../types/varinfo.hpp"

namespace tkoz
{
namespace flame
{

/*
getParams(const num_t *params, num_t& p0, num_t& p2, ...)
- store params[0],params[1],... into provided references
*/

// base case, 0 arguments
template <typename num_t>
void getParams_helper(const num_t *params, size_t n)
{
    (void)params;
    (void)n;
}

// helper with argument for array index, >=1 arguments
template <typename num_t, typename T, typename...Ts>
void getParams_helper(const num_t *params, size_t n, T& p, Ts&... ps)
{
    p = params[n]; // write current parameter
    getParams_helper(params,n+1,ps...); // increment index, write the rest
}

// put params[0],params[1],... into references (variable length)
template <typename num_t, typename...Ts>
void getParams(const num_t *params, Ts&... ps)
{
    getParams_helper(params,0,ps...);
}

/*
storeParams(std::vector<num_t>& params, p0, p1, ...)
- calls params.push_back(p) for each parameter p
*/

template <typename num_t>
void storeParams(std::vector<num_t>& params)
{
    (void)params;
}

// put parameters into the end of the parameters vector
template <typename num_t, typename T, typename...Ts>
void storeParams(std::vector<num_t>& params, T p, Ts... ps)
{
    params.push_back(p);
    storeParams(params,ps...);
}

/*
getParamsJson(const Json& j, K& k1, V& v1, K& v2, V& v2, ...)
- given JSON object, stores floating point values from keys into references
*/

void getParamsJson(const Json& j)
{
    (void)j;
}

// from JSON object, store values into references given the keys
template <typename K, typename V, typename...Ts>
void getParamsJson(const Json& j, const K& k, V& v, Ts&... ts)
{
    v = j[k].floatValue();
    getParamsJson(j,ts...);
}

// change value into  floating point type
template <typename num_t>
num_t changeType(u32 n)
{
    static_assert(std::is_same<num_t,float>::value
               || std::is_same<num_t,double>::value);
    if (std::is_same<num_t,float>::value)
    {
        union
        {
            float f;
            u32 n;
        }
        u;
        u.n = n;
        return u.f;
    }
    else
    {
        union
        {
            double f;
            u64 n;
        }
        u;
        u.n = (u64)n;
        return u.f;
    }
}

// pack 2 bytes into float/double
template <typename num_t>
num_t packBytes(u8 b0, u8 b1)
{
    u32 n0 = (u32)b0;
    u32 n1 = (u32)b1;
    return changeType<num_t>((n1 << 8) | n0);
}

// pack 3 bytes into float/double
template <typename num_t>
num_t packBytes(u8 b0, u8 b1, u8 b2)
{
    u32 n0 = (u32)b0;
    u32 n1 = (u32)b1;
    u32 n2 = (u32)b2;
    return changeType<num_t>((n2 << 16) | (n1 << 8) | n0);
}

// get 2 bytes from integer
void unpackBytes(u32 n, u8& b0, u8& b1)
{
    b0 = n & 0xFF;
    b1 = (n >>= 8) & 0xFF;
}

// get 3 bytes from integer
void unpackBytes(u32 n, u8& b0, u8& b1, u8& b2)
{
    b0 = n & 0xFF;
    b1 = (n >>= 8) & 0xFF;
    b2 = (n >>= 8) & 0xFF;
}

// get 2 bytes packed into a float
void unpackBytes(float f, u8& b0, u8& b1)
{
    union
    {
        float ff;
        u32 n;
    }
    u;
    u.ff = f;
    unpackBytes(u.n,b0,b1);
}

// get 3 bytes packed into a float
void unpackBytes(float f, u8& b0, u8& b1, u8& b2)
{
    union
    {
        float ff;
        u32 n;
    }
    u;
    u.ff = f;
    unpackBytes(u.n,b0,b1,b2);
}

// get 2 bytes packed into a double
void unpackBytes(double f, u8& b0, u8& b1)
{
    union
    {
        double ff;
        u64 n;
    }
    u;
    u.ff = f;
    unpackBytes((u32)u.n,b0,b1);
}

// get 3 bytes packed into a double
void unpackBytes(double f, u8& b0, u8& b1, u8& b2)
{
    union
    {
        double ff;
        u64 n;
    }
    u;
    u.ff = f;
    unpackBytes((u32)u.n,b0,b1,b2);
}

// convert a 2d variation function into generalized n dimensional
// the 2d variation operates on a plane specified by 2 axes
template <typename num_t, size_t dims, typename functor>
void ndFrom2dFunc(IterState<num_t,dims>& state, const num_t *params)
{
    static_assert(dims >= 2);
    if (dims == 2)
        functor(state,params);
    else // make a 2d state to call this variation function
    {
        IterState<num_t,2> state2;
        state2.rng = state.rng;
        u8 b0,b1;
        unpackBytes(params[0],b0,b1);
        state2.t = Point<num_t,2>(state.t[b0],state.t[b1]);
        state2.v = Point<num_t,2>();
        functor(state2,params+1);
        state.v[b0] += state2.v.x();
        state.v[b1] += state2.v.y();
    }
}

// convert a 3d variation function into generalized n dimensional
// the 3d variation operates on a space specified by 3 axes
template <typename num_t, size_t dims, typename functor>
void ndFrom3dFunc(IterState<num_t,dims>& state, const num_t *params)
{
    static_assert(dims >= 3);
    if (dims == 3)
        functor(state,params);
    else // make a 3d state to call this variation function
    {
        IterState<num_t,3> state3;
        state3.rng = state.rng;
        u8 b0,b1,b2;
        unpackBytes(params[0],b0,b1,b2);
        state3.t = Point<num_t,3>(state.t[b0],state.t[b1],state.t[b2]);
        state3.v = Point<num_t,3>();
        functor(state3,params+1);
        state.v[b0] += state3.v.x();
        state.v[b1] += state3.v.y();
        state.v[b2] += state3.v.z();
    }
}

// generalize a 2d parser, for higher dimensions, parse the axis parameters
template <typename num_t, size_t dims, typename functor>
void ndFrom2dParser(const Json& json, std::vector<num_t>& params)
{
    static_assert(dims >= 2);
    if (dims > 2) // axis parameters first
    {
        JsonArray axes = json["axes"];
        if (axes.size() != 2)
            throw std::runtime_error("should provide 2 axes");
        JsonInt x,y;
        x = axes[0].intValue();
        y = axes[1].intValue();
        if (x < 0 || y < 0 || x >= dims || y >= dims)
            throw std::runtime_error("axes out of range");
        params.push_back(packBytes<num_t>((u8)x,(u8)y));
    }
    functor(json,params);
}

// generalize a 3d parser, for higher dimensions, parse the axis parameters
template <typename num_t, size_t dims, typename functor>
void ndFrom3dParser(const Json& json, std::vector<num_t>& params)
{
    static_assert(dims >= 3);
    if (dims > 3) // axis parameters first
    {
        JsonArray axes = json["axes"];
        if (axes.size() != 3)
            throw std::runtime_error("should provide 3 axes");
        JsonInt x,y,z;
        x = axes[0].intValue();
        y = axes[1].intValue();
        z = axes[2].intValue();
        if (x < 0 || y < 0 || z < 0 || x >= dims || y >= dims || z >= dims)
            throw std::runtime_error("axes out of range");
        params.push_back(packBytes<num_t>((u8)x,(u8)y,(u8)z));
    }
    functor(json,params);
}

}
}
