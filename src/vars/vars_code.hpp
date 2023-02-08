/*
Variation functions
*/

#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <ctgmath>
#include <unordered_map>

#include "../types.hpp"
#include "vars_util.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz
{
namespace flame
{
namespace vars
{

/*
Linear - generalized from flam3
- use the pre-affine transformed point as is (scaled by weight)
parameters: weight
*/
template <typename num_t, size_t dims, typename rand_t>
void linear(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight = params[0];
    state.v += weight * state.t;
}

/*
Sinusoidal - generalized from flam3
- map each component to its sine
parameters: weight
*/
template <typename num_t, size_t dims, typename rand_t>
void sinusoidal(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight = params[0];
    state.v += weight * state.t.map([](num_t x){ return sin(x); });
}

/*
Spherical - generalized from flam3
- divide by squared magnitude (in 2-norm)
parameters: weight
*/
template <typename num_t, size_t dims, typename rand_t>
void spherical(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight = params[0];
    num_t r = weight / (state.t.norm2sq() + eps<num_t>::value);
    state.v += r * state.t;
}

/*
Spherical P - by TKoz, generalized further from spherical
- picks a p-norm (p > 0), spherical uses p=2
- divide by the p-th power of distance (the norm before taking the root)
parameters: weight, norm
*/
template <typename num_t, size_t dims, typename rand_t>
void spherical_p(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight,norm;
    getParams(params,weight,norm);
    num_t r = weight / (state.t.normsum(norm) + eps<num_t>::value);
    state.v += r * state.t;
}

/*
Unit Cube - by TKoz
- project point onto unit cube (scale using max norm)
parameters: weight
*/
template <typename num_t, size_t dims, typename rand_t>
void unit_cube(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight = params[0];
    num_t r = weight / (state.t.norminf() + eps<num_t>::value);
    state.v += r * state.t;
}

/*
Unit Sphere - by TKoz
- project point onto unit sphere (unit vector in same direction)
parameters: weight
*/
template <typename num_t, size_t dims, typename rand_t>
void unit_sphere(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight = params[0];
    num_t r = weight / (state.t.norm2() + eps<num_t>::value);
    state.v += r * state.t;
}

/*
Unit Sphere P - by TKoz
- project point onto unit vector in p norms
parameters: weight, norm
*/
template <typename num_t, size_t dims, typename rand_t>
void unit_sphere_p(IterState<num_t,dims,rand_t>& state, const num_t *params)
{
    num_t weight,norm;
    getParams(params,weight,norm);
    num_t r = weight / (state.t.norm(norm) + eps<num_t>::value);
    state.v += r * state.t;
}



/*

TODO REMOVE ALL THIS CODE AFTER REWRITING

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
    static void bent_parser(xform_t& xform, const Json& json,
            num_t weight, std::vector<num_t>& params)
    {
        (void)xform;
        params.push_back(weight);
        vec_t vec(json["params"]);
        for (size_t i = 0; i < dims; ++i)
            params.push_back(vec[i]);
    }
    // based on rectangles from flam3
    static void rectangles(state_t& state, params_t params)
    {
        num_t weight = params[0];
        ++params;
        vec_t vec;
        for (size_t i = 0; i < dims; ++i)
        {
            num_t p = params[i], x = state.t[i];
            if (p == 0.0)
                vec[i] = x;
            else
                vec[i] = (2.0*floor(x/p) + 1.0)*p - x;
        }
        state.v += weight * vec;
    }
    static void rectangles_parser(xform_t& xform, const Json& json,
            num_t weight, std::vector<num_t>& params)
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
    // noise from flam3 extended to higher dimensions
    static void noise(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r = weight * state.randNum();
        state.v += r * multComponents(state.t,state.randDirection());
    }
    // blur from flam3 extended to higher dimensions
    static void blur(state_t& state, params_t params)
    {
        num_t weight = params[0];
        state.v += weight * state.randNum() * state.randDirection();
    }
    // gaussian blur based on flam3 extended to higher dimensions
    static void gaussian_blur(state_t& state, params_t params)
    {
        num_t weight = params[0];
        state.v += weight * state.randGaussian() * state.randDirection();
    }
    // square from flam3
    static void square(state_t& state, params_t params)
    {
        num_t weight = params[0];
        state.v += weight * state.randPoint2();
    }
    static const data_t make_data()
    {
        data_t vardata;
        vardata["linear"] = info_t(linear);
        vardata["sinusoidal"] = info_t(sinusoidal);
        vardata["spherical"] = info_t(spherical);
        vardata["bent"] = info_t(bent,bent_parser);
        vardata["rectangles"] = info_t(rectangles,rectangles_parser);
        vardata["fisheye"] = info_t(fisheye);
        vardata["bubble"] = info_t(bubble);
        vardata["noise"] = info_t(noise);
        vardata["blur"] = info_t(blur);
        vardata["gaussian_blur"] = info_t(gaussian_blur);
        vardata["square"] = info_t(square);
        return vardata;
    }
};

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

// specialization for 2d
template <typename num_t, typename rand_t>
struct VariationsSpecific<num_t,2,rand_t>
{
    // typedefs for convenience
    typedef IterState<num_t,2,rand_t> state_t;
    typedef const num_t* params_t;
    typedef VarInfo<num_t,2,rand_t> info_t;
    typedef const XForm<num_t,2,rand_t> xform_t;
    typedef VarData<num_t,2,rand_t> data_t;
    typedef Point<num_t,2> vec_t;
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
        num_t weight = params[0] * M_1_PI; // TODO optimize, store weight/pi
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
    // from flam3, angle modified
    static void julia(state_t& state, params_t params)
    {
        num_t weight = params[0];
        static const num_t table[2] = {0.0,M_PI};
        num_t a = 0.5*state.t.angle() + table[state.randBool()];
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        state.v += state.t.norm2() * weight * vec_t(ca,sa);
    }
    // from flam3
    static void exponential(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t dx = weight * exp(state.t.x() - 1.0);
        num_t sdy,cdy;
        sincosg(M_PI*state.t.y(),&sdy,&cdy);
        state.v += dx * vec_t(cdy,sdy);
    }
    // from flam3, angle modified
    static void power(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t r,sa,ca;
        state.t.getRadiusSinCos(r,sa,ca);
        r = weight * pow(r,sa);
        state.v += r * vec_t(ca,sa);
    }
    // from flam3
    static void cosine(state_t& state, params_t params)
    {
        num_t weight = params[0];
        num_t x,y;
        state.t.getXY(x,y);
        num_t sa,ca;
        sincosg(x*M_PI,&sa,&ca);
        state.v += weight * vec_t(ca*cosh(y),-sa*sinh(y));
    }
    // from flam3, angle modified
    static void blob(state_t& state, params_t params)
    {
        num_t weight,mid,amp,waves;
        getParams(params,weight,mid,amp,waves);
        //num_t weight = params[0];
        //num_t mid = params[1];
        //num_t amp = params[2];
        //num_t waves = params[3];
        num_t r,sa,ca;
        state.t.getRadiusSinCos(r,sa,ca);
        num_t a = state.t.angle();
        r *= (mid + amp*sin(waves*a));
        state.v += weight * r * vec_t(sa,ca);
    }
    static void blob_parser(xform_t& xform, const Json& json, num_t weight,
            std::vector<num_t>& params)
    {
        (void)xform;
        num_t low,high,waves;
        getParamsJson(json,"low",low,"high",high,"waves",waves);
        // store middle point and amplitude
        storeParams(params,weight,(low+high)/2.0,(high-low)/2.0,waves);
    }
    // from flam3
    static void pdj(state_t& state, params_t params)
    {
        num_t weight,a,b,c,d;
        getParams(params,weight,a,b,c,d);
        num_t nx1 = cos(b*state.t.x());
        num_t nx2 = sin(c*state.t.x());
        num_t ny1 = sin(a*state.t.y());
        num_t ny2 = cos(d*state.t.y());
        state.v += weight * vec_t(ny1-nx1,nx2-ny2);
    }
    static void pdj_parser(xform_t& xform, const Json& json, num_t weight,
            std::vector<num_t>& params)
    {
        (void)xform;
        num_t a,b,c,d;
        getParamsJson(json,"a",a,"b",b,"c",c,"d",d);
        storeParams(params,weight,a,b,c,d);
    }
    // from flam3
    static void cylinder(state_t& state, params_t params)
    {
        num_t weight = params[0];
        state.v += weight * vec_t(sin(state.t.x()),state.t.y());
    }
    // from flam3
    static void perspective(state_t& state, params_t params)
    {
        num_t weight,dist,vsin,vfcos;
        getParams(params,weight,dist,vsin,vfcos);
        num_t wt = weight / (dist - state.t.y()*vsin);
        state.v += wt * vec_t(dist*state.t.x(),vfcos*state.t.y());
    }
    static void perspective_parser(xform_t& xform, const Json& json,
            num_t weight, std::vector<num_t>& params)
    {
        (void)xform;
        num_t dist,angle;
        getParamsJson(json,"distance",dist,"angle",angle);
        storeParams(params,weight,dist,sin(angle),dist*cos(angle));
    }
    // from flam3
    static void julian(state_t& state, params_t params)
    {
        num_t weight,abspower,invpower,cn;
        getParams(params,weight,abspower,invpower,cn);
        i32 t = trunc(abspower*state.randNum());
        num_t a = (state.t.angle() + (2.0*M_PI)*t) * invpower;
        num_t r = weight * pow(state.t.norm2sq(),cn);
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        state.v += r * vec_t(ca,sa);
    }
    static void julian_parser(xform_t& xform, const Json& json,
            num_t weight, std::vector<num_t>& params)
    {
        (void)xform;
        num_t power,dist;
        getParamsJson(json,"power",power,"distance",dist);
        storeParams(params,weight,fabs(power),1.0/power,dist/(2.0*power));
    }
    // from flam3
    static void juliascope(state_t& state, params_t params)
    {
        num_t weight,abspower,invpower,cn;
        getParams(params,weight,abspower,invpower,cn);
        i32 t = trunc(abspower*state.randNum());
        static const num_t sign[2] = {-1.0,1.0};
        num_t a = ((2.0*M_PI)*t + sign[t&1]*state.t.angle()) * invpower;
        num_t r = weight * pow(state.t.norm2sq(),cn);
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        state.v += r * vec_t(ca,sa);
    }
    static void juliascope_parser(xform_t& xform, const Json& json,
            num_t weight, std::vector<num_t>& params)
    {
        (void)xform;
        num_t power,dist;
        getParamsJson(json,"power",power,"distance",dist);
        storeParams(params,weight,fabs(power),1.0/power,dist/(2.0*power));
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
        vardata["julia"] = info_t(julia);
        vardata["exponential"] = info_t(exponential);
        vardata["power"] = info_t(power);
        vardata["cosine"] = info_t(cosine);
        vardata["blob"] = info_t(blob,blob_parser);
        vardata["pdj"] = info_t(pdj,pdj_parser);
        vardata["cylinder"] = info_t(cylinder);
        vardata["perspective"] = info_t(perspective,perspective_parser);
        vardata["julian"] = info_t(julian,julian_parser);
        vardata["juliascope"] = info_t(juliascope,juliascope_parser);
        return vardata;
    }
};

*/

}
}
}

#undef likely
#undef unlikely
