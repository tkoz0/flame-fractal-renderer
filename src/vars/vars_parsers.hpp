/*
Functions for parsing variation parameters
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

#include "../types.hpp"
#include "vars_util.hpp"

namespace tkoz
{
namespace flame
{
namespace parsers
{

/*
Weight only - for variations using only 1 weight parameter
*/
template <typename num_t>
static void weight_only(const Json& json, std::vector<num_t>& params)
{
    params.push_back(json["weight"].floatValue());
}

/*
Spherical P - parse weight and the norm parameter
*/
template <typename num_t>
static void spherical_p(const Json& json, std::vector<num_t>& params)
{
    weight_only(json,params);
    num_t p = json["norm"].floatValue();
    if (p <= 0.0)
        throw std::runtime_error("spherical_p norm not positive");
    params.push_back(p);
}

}
}
}
