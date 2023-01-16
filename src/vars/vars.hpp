
#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <ctgmath>
#include <unordered_map>

#include "../types.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define ENABLE_IF(COND,RET) template <typename RET2 = RET> \
    typename std::enable_if<(COND),RET2>::type
#define ENABLE_IFEQ(N1,N2,RET) template <typename RET2 = RET> \
    typename tkoz::flame::enable_if_eq<N1,N2,RET2>::type

namespace tkoz
{
namespace flame
{

// define getParams function for getting numbers from params pointer

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

// push_back multiple values to the params vector
template <typename num_t>
void storeParams(std::vector<num_t>& params)
{
    (void)params;
}

template <typename num_t, typename T, typename...Ts>
void storeParams(std::vector<num_t>& params, T p, Ts... ps)
{
    params.push_back(p);
    storeParams(params,ps...);
}

// extract values from json to store in referencesk
void getParamsJson(const Json& j)
{
    (void)j;
}

template <typename K, typename V, typename...Ts>
void getParamsJson(const Json& j, const K& k, V& v, Ts&... ts)
{
    v = j[k].floatValue();
    getParamsJson(j,ts...);
}

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
struct Variations
{
    static const VarData<num_t,dims,rand_t> vars;
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
struct VariationsGeneric
{
    // typedefs for convenience
    typedef IterState<num_t,dims,rand_t> state_t;
    typedef const num_t* params_t;
    typedef VarInfo<num_t,dims,rand_t> info_t;
    typedef const XForm<num_t,dims,rand_t> xform_t;
    typedef VarData<num_t,dims,rand_t> data_t;
    typedef Point<num_t,dims> vec_t;
    static constexpr num_t eps = tkoz::flame::eps<num_t>::value;
    // from flam3
    static void linear(state_t& state, params_t params)
    {
        num_t weight = params[0];
        state.v += weight * state.t;
    }
    // from flam3
    static void sinusoidal(state_t& state, params_t params)
    {
        num_t weight = params[0];
        state.v += weight * state.t.map([](num_t x){ return sin(x); });
    }
    // from flam3
    static void spherical(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = weight / (state.t.norm2sq() + eps);
        state.v += r * state.t;
    }
    // based on bent2 from flam3
    static void bent(state_t& state, params_t params)
    {
        num_t weight = params[0];
        ++params;
        vec_t vec = state.t;
        for (size_t i = 0; i < dims; ++i)
            if (vec[i] < 0)
                vec[i] *= params[i];
        state.v += weight * vec;
    }
    static void bent_parser(xform_t& xform, const Json& json, num_t weight,
            std::vector<num_t>& params)
    {
        (void)xform;
        params.push_back(weight);
        vec_t vec(json["params"]);
        for (size_t i = 0; i < dims; ++i)
            params.push_back(vec[i]);
    }
    // fisheye based on flam3 with order corrected, constant changed
    static void fisheye(state_t& state, params_t params)
    {
        num_t weight = params[0];
        // TODO support changing the 1.0 constant
        num_t r = weight / (state.t.norm2() + 1.0);
        state.v += r * state.t;
    }
    // bubble based on flam3, constant changed
    static void bubble(state_t& state, params_t params)
    {
        num_t weight = params[0];
        // TODO support changing the 4.0 constant
        num_t r = weight / (state.t.norm2sq() + 4.0);
        state.v += r * state.t;
    }
    // noise from flam3
    ENABLE_IFEQ(dims,2,void) static noise(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = weight * state.randNum();
        num_t a = (2.0*M_PI) * state.randNum();
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        state.v += r * multComponents(state.t,vec_t(ca,sa));
    }
    // noise extended to 3d
    ENABLE_IFEQ(dims,3,void) noise(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = weight * state.randNum();
        num_t p = acos(2.0*state.randNum()-1.0);
        num_t t = (2.0*M_PI) * state.randNum();
        // coordinates on unit sphere
        num_t x = sin(t)*cos(p);
        num_t y = sin(t)*sin(p);
        num_t z = cos(t);
        state.v += r * multComponents(state.t,vec_t(x,y,z));
    }
    // noise extended to higher dimensions
    ENABLE_IF(dims>3,void) static noise(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = weight * state.randNum();
        Point<num_t,dims> a = state.randAngle();
        state.v += r * multComponents(state.t,a);
    }
    static const data_t make_data()
    {
        data_t vardata;
        vardata["linear"] = info_t(linear);
        vardata["sinusoidal"] = info_t(sinusoidal);
        vardata["spherical"] = info_t(spherical);
        vardata["bent"] = info_t(bent,bent_parser);
        vardata["fisheye"] = info_t(fisheye);
        vardata["bubble"] = info_t(bubble);
        return vardata;
    }
};

// dimension specific variations
// specialize this class in separate files per dimension
template <typename num_t, size_t dims, typename rand_t>
struct VariationsSpecific;

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

#undef ENABLE_IF
#undef ENABLE_IFEQ
