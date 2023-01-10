
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

// specialization for 2d
template <typename num_t, typename rand_t>
class VariationsSpecific<num_t,2,rand_t>
{
private:
    // from flam3
    static void swirl(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t sr,cr;
        sincosg(state.t.norm2sq(),&sr,&cr);
        num_t x = state.t.x(), y = state.t.y();
        state.v += weight * Point<num_t,2>(sr*x-cr*y,cr*x+sr*y);
    }
    // from flam3
    static void horseshoe(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t r = weight / (state.t.norm2() + eps<num_t>::value);
        num_t x = state.t.x(), y = state.t.y();
        state.v += r * Point<num_t,2>((x-y)*(x+y),2.0*x*y);
    }
    // from flam3, angle modified
    static void polar(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = atan2(state.t.y(),state.t.x());
        state.v += weight * Point<num_t,2>(a*M_1_PI,state.t.norm2()-1.0);
    }
    // from flam3, angle modified
    static void handkerchief(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = atan2(state.t.y(),state.t.x());
        num_t r = state.t.norm2();
        state.v += weight * r * Point<num_t,2>(sin(a+r),cos(a-r));
    }
    // from flam3, angle modified
    static void heart(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t r = state.t.norm2();
        num_t a = atan2(state.t.y(),state.t.x());
        num_t sa,ca;
        sincosg(r*a,&sa,&ca);
        state.v += weight * r * Point<num_t,2>(sa,-ca);
    }
    // from flam3, angle modified
    static void disc(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0] * M_1_PI; // TODO optimize
        num_t a = atan2(state.t.y(),state.t.x()) * weight;
        num_t r = state.t.norm2();
        num_t sr,cr;
        sincosg(M_PI*r,&sr,&cr);
        state.v += a * Point<num_t,2>(sr,cr);
    }
    // from flam3, angle modified
    static void spiral(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = atan2(state.t.y(),state.t.x());
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        num_t r = state.t.norm2() + eps<num_t>::value;
        num_t sr,cr;
        sincosg(r,&sr,&cr);
        num_t r1 = weight / r;
        state.v += r1 * Point<num_t,2>(ca+sr,sa-cr);
    }
public:
    static const VarData<num_t,2,rand_t> make_data()
    {
        VarData<num_t,2,rand_t> vardata;
        vardata["swirl"] = VarInfo<num_t,2,rand_t>(swirl);
        vardata["horseshoe"] = VarInfo<num_t,2,rand_t>(horseshoe);
        vardata["polar"] = VarInfo<num_t,2,rand_t>(polar);
        vardata["handkerchief"] = VarInfo<num_t,2,rand_t>(handkerchief);
        vardata["heart"] = VarInfo<num_t,2,rand_t>(heart);
        vardata["disc"] = VarInfo<num_t,2,rand_t>(disc);
        vardata["spiral"] = VarInfo<num_t,2,rand_t>(spiral);
        return vardata;
    }
};

}
}

#undef likely
#undef unlikely
