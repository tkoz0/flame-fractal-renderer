/*
Helper functions for dealing with variations
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

#include "../types/types.hpp"

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

}
}
