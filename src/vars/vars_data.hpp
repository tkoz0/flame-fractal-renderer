
#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>
#include <unordered_map>

#include "../types.hpp"
#include "vars_code.hpp"
#include "vars_parsers.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz
{
namespace flame
{

template <typename num_t, size_t dims, typename rand_t>
using VarFunc = std::function<void(IterState<num_t,dims,rand_t>&,
                              const num_t*)>;

template <typename num_t, size_t dims, typename rand_t>
using VarParse = std::function<void(const Json&,std::vector<num_t>&)>;

// representation of variations data
template <typename num_t, size_t dims, typename rand_t>
class VarInfo
{
private:
    // function pointer for calculating the variation
    VarFunc<num_t,dims,rand_t> func;
    // params parser, use nullptr for default of weight only
    VarParse<num_t,dims,rand_t> params;
public:
    // default constructable for use with unordered_map operator[]
    VarInfo(VarFunc<num_t,dims,rand_t> func = nullptr,
            VarParse<num_t,dims,rand_t> params = nullptr):
        func(func),params(params) {}
    inline VarFunc<num_t,dims,rand_t> getFPtr() const
    {
        return func;
    }
    inline VarParse<num_t,dims,rand_t> getPPtr() const
    {
        return params;
    }
};

template <typename num_t, size_t dims, typename rand_t>
using VarData = std::unordered_map<std::string,VarInfo<num_t,dims,rand_t>>;

template <typename num_t, size_t dims, typename rand_t>
struct Variations
{
    // typedefs for convenience
    typedef IterState<num_t,dims,rand_t> state_t;
    typedef const num_t* params_t;
    typedef VarInfo<num_t,dims,rand_t> info_t;
    typedef const XForm<num_t,dims,rand_t> xform_t;
    typedef VarData<num_t,dims,rand_t> data_t;
    typedef Point<num_t,dims> vec_t;
    static constexpr num_t eps = eps<num_t>::value;
    // variations data
    static const VarData<num_t,dims,rand_t> vars;
    // setup the variations data
    static VarData<num_t,dims,rand_t> make_map()
    {
        VarData<num_t,dims,rand_t> data;
        data["linear"] = info_t(
            vars::linear<num_t,dims,rand_t>,
            parsers::weight_only<num_t>
        );
        data["sinusoidal"] = info_t(
            vars::sinusoidal<num_t,dims,rand_t>,
            parsers::weight_only<num_t>
        );
        data["spherical"] = info_t(
            vars::spherical<num_t,dims,rand_t>,
            parsers::weight_only<num_t>
        );
        data["spherical_p"] = info_t(
            vars::spherical_p<num_t,dims,rand_t>,
            parsers::spherical_p<num_t>
        );
        data["unit_cube"] = info_t(
            vars::unit_cube<num_t,dims,rand_t>,
            parsers::weight_only<num_t>
        );
        data["unit_sphere"] = info_t(
            vars::unit_sphere<num_t,dims,rand_t>,
            parsers::weight_only<num_t>
        );
        data["unit_sphere_p"] = info_t(
            vars::unit_sphere_p<num_t,dims,rand_t>,
            parsers::spherical_p<num_t>
        );
        return data;
    }
    // access variations data based on its name
    static const VarInfo<num_t,dims,rand_t>& get(const std::string& s)
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
template <typename num_t, size_t dims, typename rand_t>
const VarData<num_t,dims,rand_t> Variations<num_t,dims,rand_t>::vars =
    Variations<num_t,dims,rand_t>::make_map();

}
}

#undef likely
#undef unlikely
