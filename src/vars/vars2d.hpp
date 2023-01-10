
#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <ctgmath>
#include <unordered_map>

#include "vars.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz
{
namespace flame
{

template <typename num_t, typename rand_t>
class Variations<num_t,2,rand_t>
{
private:
    static const VarData<num_t,2,rand_t> vars;
public:
    static const VarInfo<num_t,2,rand_t>& get(const std::string& s)
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
};

template <typename num_t, typename rand_t>
struct vars2d { static const VarData<num_t,2,rand_t> data; };

template <typename num_t, typename rand_t>
const VarData<num_t,2,rand_t> vars2d<num_t,rand_t>::data = {
    {"swirl",VarInfo<num_t,2,rand_t>(
        [](IterState<num_t,2,rand_t>& state, const num_t *params)
        {
            num_t weight = params[0];
            num_t sr,cr;
            sincosg(state.t.norm2sq(),&sr,&cr);
            num_t x = state.t.x(), y = state.t.y();
            state.v += weight * Point<num_t,2>(sr*x-cr*y,cr*x+sr*y);
        }
    )},
    {"horseshoe",VarInfo<num_t,2,rand_t>(
        [](IterState<num_t,2,rand_t>& state, const num_t *params)
        {
            num_t weight = params[0];
            num_t r = weight / (state.t.norm2() + eps<num_t>::value);
            num_t x = state.t.x(), y = state.t.y();
            state.v += r * Point<num_t,2>((x-y)*(x+y),2.0*x*y);
        }
    )}
};

template <typename num_t, typename rand_t>
VarData<num_t,2,rand_t> make_vars2d()
{
    VarData<num_t,2,rand_t> data;
    add_generic_variations(data);
    for (auto itr : vars2d<num_t,rand_t>::data)
        data.insert(itr);
    return data;
}

template <typename num_t, typename rand_t>
const VarData<num_t,2,rand_t> Variations<num_t,2,rand_t>::vars =
    make_vars2d<num_t,rand_t>();

}
}

#undef N
#undef likely
#undef unlikely
