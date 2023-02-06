/*
Functions for parsing variation parameters
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>

#include "../types.hpp"

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

}
}
}