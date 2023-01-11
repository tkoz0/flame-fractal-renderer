
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
struct VariationsSpecific<num_t,2,rand_t>
{
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
        num_t a = state.t.angle();
        state.v += weight * Point<num_t,2>(a*M_1_PI,state.t.norm2()-1.0);
    }
    // from flam3, angle modified
    static void handkerchief(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        num_t r = state.t.norm2();
        state.v += weight * r * Point<num_t,2>(sin(a+r),cos(a-r));
    }
    // from flam3, angle modified
    static void heart(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t r = state.t.norm2();
        num_t a = state.t.angle();
        num_t sa,ca;
        sincosg(r*a,&sa,&ca);
        state.v += weight * r * Point<num_t,2>(sa,-ca);
    }
    // from flam3, angle modified
    static void disc(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0] * M_1_PI; // TODO optimize
        num_t a = state.t.angle() * weight;
        num_t r = state.t.norm2();
        num_t sr,cr;
        sincosg(M_PI*r,&sr,&cr);
        state.v += a * Point<num_t,2>(sr,cr);
    }
    // from flam3, angle modified
    static void spiral(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        num_t r = state.t.norm2() + eps<num_t>::value;
        num_t sr,cr;
        sincosg(r,&sr,&cr);
        num_t r1 = weight / r;
        state.v += r1 * Point<num_t,2>(ca+sr,sa-cr);
    }
    // from flam3, angle modified
    static void hyperbolic(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t r = state.t.norm2() + eps<num_t>::value;
        num_t a = state.t.angle();
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        state.v += weight * Point<num_t,2>(sa/r,ca*r);
    }
    // from flam3, angle modified
    static void diamond(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        num_t r = state.t.norm2();
        num_t sr,cr;
        sincosg(r,&sr,&cr);
        state.v += weight * Point<num_t,2>(sa*cr,ca*sr);
    }
    // from flam3, angle modified
    static void ex(IterState<num_t,2,rand_t>& state, const num_t *params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        num_t r = state.t.norm2();
        num_t n0 = sin(a+r);
        num_t n1 = cos(a-r);
        num_t m0 = n0*n0*n0 * r;
        num_t m1 = n1*n1*n1 * r;
        state.v += weight * Point<num_t,2>(m0+m1,m0-m1);
    }
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
        vardata["hyperbolic"] = VarInfo<num_t,2,rand_t>(hyperbolic);
        vardata["diamond"] = VarInfo<num_t,2,rand_t>(diamond);
        vardata["ex"] = VarInfo<num_t,2,rand_t>(ex);
        return vardata;
    }
};

}
}

#undef likely
#undef unlikely
