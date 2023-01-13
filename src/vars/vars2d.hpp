
#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <ctgmath>
#include <unordered_map>

#include "vars.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)
#define DIMS 2

namespace tkoz
{
namespace flame
{

// specialization for 2d
template <typename num_t, typename rand_t>
struct VariationsSpecific<num_t,DIMS,rand_t>
{
    // typedefs for convenience
    typedef IterState<num_t,DIMS,rand_t> state_t;
    typedef const num_t* params_t;
    typedef VarInfo<num_t,DIMS,rand_t> info_t;
    typedef const XForm<num_t,DIMS,rand_t> xform_t;
    typedef VarData<num_t,DIMS,rand_t> data_t;
    typedef Point<num_t,DIMS> vec_t;
    static constexpr num_t eps = tkoz::flame::eps<num_t>::value;
    // from flam3
    static void swirl(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t sr,cr;
        sincosg(state.t.norm2sq(),&sr,&cr);
        num_t x,y;
        state.t.getXY(x,y);
        state.v += weight * vec_t(sr*x-cr*y,cr*x+sr*y);
    }
    // from flam3
    static void horseshoe(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = weight / (state.t.norm2() + eps);
        num_t x,y;
        state.t.getXY(x,y);
        state.v += r * vec_t((x-y)*(x+y),2.0*x*y);
    }
    // from flam3, angle modified
    static void polar(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        state.v += weight * vec_t(a*M_1_PI,state.t.norm2()-1.0);
    }
    // from flam3, angle modified
    static void handkerchief(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        num_t r = state.t.norm2();
        state.v += weight * r * vec_t(sin(a+r),cos(a-r));
    }
    // from flam3, angle modified
    static void heart(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = state.t.norm2();
        num_t a = state.t.angle();
        num_t sa,ca;
        sincosg(r*a,&sa,&ca);
        state.v += weight * r * vec_t(sa,-ca);
    }
    // from flam3, angle modified
    static void disc(state_t& state, params_t params)
    {
        num_t weight = params[0] * M_1_PI; // TODO optimize
        num_t a = state.t.angle() * weight;
        num_t r = state.t.norm2();
        num_t sr,cr;
        sincosg(M_PI*r,&sr,&cr);
        state.v += a * vec_t(sr,cr);
    }
    // from flam3, angle modified
    static void spiral(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r,sa,ca;
        state.t.getRadiusSinCos(r,sa,ca);
        r += eps;
        num_t sr,cr;
        sincosg(r,&sr,&cr);
        num_t r1 = weight / r;
        state.v += r1 * vec_t(ca+sr,sa-cr);
    }
    // from flam3, angle modified
    static void hyperbolic(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r,sa,ca;
        state.t.getRadiusSinCos(r,sa,ca);
        r += eps;
        state.v += weight * vec_t(sa/r,ca*r);
    }
    // from flam3, angle modified
    static void diamond(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r,sa,ca;
        state.t.getRadiusSinCos(r,sa,ca);
        num_t sr,cr;
        sincosg(r,&sr,&cr);
        state.v += weight * vec_t(sa*cr,ca*sr);
    }
    // from flam3, angle modified
    static void ex(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t a = state.t.angle();
        num_t r = state.t.norm2();
        num_t n0 = sin(a+r);
        num_t n1 = cos(a-r);
        num_t m0 = n0*n0*n0 * r;
        num_t m1 = n1*n1*n1 * r;
        state.v += weight * vec_t(m0+m1,m0-m1);
    }
    static const data_t make_data()
    {
        data_t vardata;
        vardata["swirl"] = info_t(swirl);
        vardata["horseshoe"] = info_t(horseshoe);
        vardata["polar"] = info_t(polar);
        vardata["handkerchief"] = info_t(handkerchief);
        vardata["heart"] = info_t(heart);
        vardata["disc"] = info_t(disc);
        vardata["spiral"] = info_t(spiral);
        vardata["hyperbolic"] = info_t(hyperbolic);
        vardata["diamond"] = info_t(diamond);
        vardata["ex"] = info_t(ex);
        return vardata;
    }
};

}
}

#undef likely
#undef unlikely
#undef DIMS
