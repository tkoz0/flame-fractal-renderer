
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
class VarInfo
{
private:
    // function pointer for calculating the variation
    VarFunc<num_t,dims,rand_t> func;
    // params parser, use nullptr for default of weight only
    VarParse<num_t,dims,rand_t> params;
    // precalculate flags (not used currently)
    u32 pc_flags;
public:
    // default constructable for use with unordered_map operator[]
    VarInfo(VarFunc<num_t,dims,rand_t> func = nullptr,
            VarParse<num_t,dims,rand_t> params = nullptr,
            u32 pc_flags = 0):
        func(func),params(params),pc_flags(pc_flags) {}
    inline VarFunc<num_t,dims,rand_t> getFPtr() const { return func; }
    inline VarParse<num_t,dims,rand_t> getPPtr() const { return params; }
    inline u32 getPCFlags() const { return pc_flags; }
};

template <typename num_t, size_t dims, typename rand_t>
using VarData = std::unordered_map<std::string,VarInfo<num_t,dims,rand_t>>;

template <typename num_t, size_t dims, typename rand_t>
class Variations
{
private:
    static const VarData<num_t,dims,rand_t> vars;
public:
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
    static size_t size()
    {
        return vars.size();
    }
};

// dimension generic variations
template <typename num_t, size_t dims, typename rand_t>
class VariationsGeneric
{
private:
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
public:
    static const VarData<num_t,dims,rand_t> make_data()
    {
        VarData<num_t,dims,rand_t> vardata;
        vardata["linear"] = VarInfo<num_t,dims,rand_t>(linear);
        vardata["sinusoidal"] = VarInfo<num_t,dims,rand_t>(sinusoidal);
        vardata["spherical"] = VarInfo<num_t,dims,rand_t>(spherical);
        return vardata;
    }
};

// dimension specific variations
// specialize this class in separate files per dimension
template <typename num_t, size_t dims, typename rand_t>
class VariationsSpecific;

// combine dimension generic and dimension specific variations
template <typename K, typename V, typename ...Maps>
std::unordered_map<K,V> combine_maps(const Maps&... maps)
{
    std::unordered_map<K,V> ret;
    for (auto map : {maps...})
        for (auto itr : map)
            ret.insert(itr);
    return ret;
}

template <typename num_t, size_t dims, typename rand_t>
const VarData<num_t,dims,rand_t> Variations<num_t,dims,rand_t>::vars =
    combine_maps<std::string,VarInfo<num_t,dims,rand_t>,
            VarData<num_t,dims,rand_t>>(
        VariationsGeneric<num_t,dims,rand_t>::make_data(),
        VariationsSpecific<num_t,dims,rand_t>::make_data());

}
}

#undef likely
#undef unlikely
