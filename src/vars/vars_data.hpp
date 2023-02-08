/*
Class for hardcoded variations information
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>
#include <unordered_map>

#include "../types/types.hpp"
#include "../types/xform.hpp"
#include "../types/varinfo.hpp"
#include "vars_code.hpp"
#include "vars_parsers.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz
{
namespace flame
{

template <typename num_t, size_t dims>
struct Variations
{
    // typedefs for convenience
    typedef IterState<num_t,dims> state_t;
    typedef const num_t* params_t;
    typedef VarInfo<num_t,dims> info_t;
    typedef const XForm<num_t,dims> xform_t;
    typedef VarData<num_t,dims> data_t;
    typedef Point<num_t,dims> vec_t;
    static constexpr num_t eps = eps<num_t>::value;
    // variations data
    static const VarData<num_t,dims> vars;
    // setup the variations data
    static VarData<num_t,dims> make_map()
    {
        VarData<num_t,dims> data;
        data["linear"] = info_t(
            vars::linear<num_t,dims>,
            parsers::weight_only<num_t>
        );
        data["sinusoidal"] = info_t(
            vars::sinusoidal<num_t,dims>,
            parsers::weight_only<num_t>
        );
        data["spherical"] = info_t(
            vars::spherical<num_t,dims>,
            parsers::weight_only<num_t>
        );
        data["spherical_p"] = info_t(
            vars::spherical_p<num_t,dims>,
            parsers::spherical_p<num_t>
        );
        data["unit_cube"] = info_t(
            vars::unit_cube<num_t,dims>,
            parsers::weight_only<num_t>
        );
        data["unit_sphere"] = info_t(
            vars::unit_sphere<num_t,dims>,
            parsers::weight_only<num_t>
        );
        data["unit_sphere_p"] = info_t(
            vars::unit_sphere_p<num_t,dims>,
            parsers::spherical_p<num_t>
        );
        return data;
    }
    // access variations data based on its name
    static const VarInfo<num_t,dims>& get(const std::string& s)
    {
        try
        {
            return vars.at(s);
        }
        catch (std::out_of_range& ex)
        {
            throw std::runtime_error("invalid variation: "+s);
        }
    }
    // number of variations
    static size_t size()
    {
        return vars.size();
    }
};

// assign the hardcoded data map to the static variable
template <typename num_t, size_t dims>
const VarData<num_t,dims> Variations<num_t,dims>::vars =
    Variations<num_t,dims>::make_map();

}
}

#undef likely
#undef unlikely
