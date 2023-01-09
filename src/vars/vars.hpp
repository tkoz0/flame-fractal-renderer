
#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <ctgmath>
#include <unordered_map>

#include "../types.hpp"

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
using VarParse = std::function<void(const XForm<num_t,dims,rand_t>&,
                                    const Json&,
                                    num_t,
                                    std::vector<num_t>&)>;

// representation of variations data
template <typename num_t, size_t dims, typename rand_t>
struct VarInfo
{
    // function pointer for calculating the variation
    const VarFunc<num_t,dims,rand_t> func;
    // precalculate flags
    const u32 pc_flags;
    // params parser, use nullptr for default of weight only
    const VarParse<num_t,dims,rand_t> params;
    VarInfo(VarFunc<num_t,dims,rand_t> func,
            VarParse<num_t,dims,rand_t> params = nullptr,
            u32 pc_flags = 0):
        func(func),pc_flags(pc_flags),params(params) {}
};

template <typename num_t, size_t dims, typename rand_t>
using VarData = std::unordered_map<std::string,VarInfo<num_t,dims,rand_t>>;

template <typename num_t, size_t dims, typename rand_t>
class Variations;

// dimension generic variations
template <typename num_t, size_t dims, typename rand_t>
class VariationsGeneric
{
public:
    static void linear(IterState<num_t,dims,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        state.v += weight * state.t;
    }
    static void sinusoidal(IterState<num_t,dims,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        state.v += weight * state.t.map([](num_t x){ return sin(x); });
    }
    static void spherical(IterState<num_t,dims,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t r = weight / (state.t.norm2sq() + eps<num_t>::value);
        state.v += r * state.t;
    }
};

// add generic variations to classes for each number of dimensions
template <typename num_t, size_t dims, typename rand_t>
void add_generic_variations(VarData<num_t,dims,rand_t>& vardata)
{
    vardata.insert(std::make_pair("linear",VarInfo<num_t,dims,rand_t>(
        VariationsGeneric<num_t,dims,rand_t>::linear
    )));
    vardata.insert(std::make_pair("sinusoidal",VarInfo<num_t,dims,rand_t>(
        VariationsGeneric<num_t,dims,rand_t>::sinusoidal
    )));
    vardata.insert(std::make_pair("spherical",VarInfo<num_t,dims,rand_t>(
        VariationsGeneric<num_t,dims,rand_t>::spherical
    )));
}

}
}

#undef likely
#undef unlikely
