#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <ctgmath>
#include <unordered_map>

#include "types.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

// types
#define VAR_T   VarInfo<num_t,rand_t>
#define STATE_T IterState<num_t,rand_t>&
#define XFORM_T const XForm<num_t,rand_t>&
#define JSON_T  const Json&
#define PARAM_T std::vector<num_t>&
#define VEC_T   Point2D<num_t>

// macros for variation functions
//#define VAR_FUNC  [](STATE_T S, num_t W, const num_t *P)
#define VAR_FUNC  [](STATE_T state, const num_t *params)
#define VAR_PARSE [](XFORM_T xform, JSON_T json, num_t weight, PARAM_T varp)
#define VAR_RET(ret) state.v += (ret)
#define VEC(x,y) VEC_T(x,y)
#define TX state.t.x
#define TY state.t.y
#define TP state.t
#define EPS eps<num_t>::value

namespace tkoz
{
namespace flame
{

// iterator state used by variation functions, defined with renderer
template <typename num_t, typename rand_t> struct IterState;

// variation info for parsing variations
template <typename num_t, typename rand_t> struct VarInfo
{
    const std::function<void(STATE_T,const num_t*)> func;
    const std::function<void(XFORM_T,JSON_T,num_t,PARAM_T)> params;
    VarInfo(std::function<void(STATE_T,const num_t*)> func,
        std::function<void(XFORM_T,JSON_T,num_t,PARAM_T)> params = nullptr):
        func(func),params(params) {}
};

// struct containing hardcoded information of available variations
template <typename num_t, typename rand_t> struct vars
{
    static const std::unordered_map<std::string,VAR_T> data;
};

// parse JSON or use default value
template <typename num_t>
num_t parse_var_param(const Json& j, const std::string& key,
        num_t default_, bool allow_integer = false)
{
    Json value;
    if (j.valueAt(key,value))
    {
        if (value.isFloat()) return value.floatValue();
        if (value.isInt()) return (num_t) value.intValue();
        if (value.isNull()) return default_;
        throw std::runtime_error("value is not a number");
    }
    else
        return default_;
}

/*
Variation functions, parameters IterState<num_t>& state, const num_t *params
- Inputs are the iteration state and pointer to parameters
- Each computes a point (vector) from S.t to add to S.v (the variation sum)
  - use the TX,TY,TP,VAR_RET macros above
- Additional details and the RNG are accessible through state
- Parameters and some precomputed values are accessible through params
- if the parser function is not specified
  - params[0] is the weight (0 should have the effect of S.v += (0,0))

Variations using extra parameters have another function to create them
- They call push_back to the vector `varp` for each precomputed parameter
- Parameters can be parsed from the JSON data for the variation
- Details about the xform can be used (such as the affine transforms)
- If not specified, then params[0] is weight and there are no others

TODO precomputed variables based on S.t for each iteration
*/

template <typename num_t, typename rand_t>
const std::unordered_map<std::string,VAR_T>
vars<num_t,rand_t>::data =
{
    {"linear", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            VAR_RET(W * TP);
        }
    )},
    {"sinusoidal", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            VAR_RET(W * VEC(sin(TX),sin(TY)));
        }
    )},
    {"spherical", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r = W / (TP.r2() + EPS);
            VAR_RET(r * TP);
        }
    )},
    {"swirl", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t sr,cr;
            sincosg(TP.r2(),&sr,&cr);
            VAR_RET(W * VEC(sr*TX-cr*TY,cr*TX+sr*TY));
        }
    )},
    {"horseshoe", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r = W / (TP.r() + EPS);
            VAR_RET(r * VEC((TX-TY)*(TX+TY),2.0*TX*TY));
        }
    )},
    {"polar", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            VAR_RET(W * VEC(TP.atan()*M_1_PI,TP.r()-1.0));
        }
    )},
    {"handkerchief", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = TP.atan();
            num_t r = TP.r();
            VAR_RET(W * r * VEC(sin(a+r),cos(a-r)));
        }
    )},
    {"heart", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r = TP.r();
            num_t sa,ca;
            sincosg(r*TP.atan(),&sa,&ca);
            VAR_RET(r * W * VEC(sa,-ca));
        }
    )},
    {"disc", VAR_T(
        VAR_FUNC
        {
            num_t W_pi = params[0]; // W/pi
            num_t a = TP.atan() * W_pi;
            num_t sr,cr;
            sincosg(M_PI*TP.r(),&sr,&cr);
            VAR_RET(a * VEC(sr,cr));
        },
        // store: weight/pi
        VAR_PARSE
        {
            varp.push_back(M_1_PI*weight);
        }
    )},
    {"spiral", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t sa,ca;
            sincosg(TP.atan(),&sa,&ca);
            num_t r = TP.r() + EPS;
            num_t sr,cr;
            sincosg(r,&sr,&cr);
            num_t r1 = W / r;
            VAR_RET(r1 * VEC(ca+sr,sa-cr));
        }
    )},
    {"hyperbolic", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r = TP.r() + EPS;
            num_t sa,ca;
            sincosg(TP.atan(),&sa,&ca);
            VAR_RET(W * VEC(sa/r,ca*r));
        }
    )},
    {"diamond", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t sa,ca;
            sincosg(TP.atan(),&sa,&ca);
            num_t sr,cr;
            sincosg(TP.r(),&sr,&cr);
            VAR_RET(W * VEC(sa*cr,ca*sr));
        }
    )},
    {"ex", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = TP.atan();
            num_t r = TP.r();
            num_t n0 = sin(a+r);
            num_t n1 = cos(a-r);
            // TODO is this the best way to compute cubes
            num_t m0 = n0*n0*n0 * r;
            num_t m1 = n1*n1*n1 * r;
            VAR_RET(W * VEC(m0+m1,m0-m1));
        }
    )},
    {"julia", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            static const num_t table[2] = {0.0,M_PI};
            num_t r = TP.r() * W;
            num_t sa,ca;
            sincosg(0.5*TP.atan()+table[state.randBool()],&sa,&ca);
            VAR_RET(r * VEC(ca,sa));
        }
    )},
    {"bent", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            static const num_t table_x[2] = {1.0,2.0};
            static const num_t table_y[2] = {1.0,0.5};
            num_t x = TX;
            num_t y = TY;
            x *= table_x[x < 0.0];
            y *= table_y[y < 0.0];
            VAR_RET(W * VEC(x,y));
        }
    )},
    {"waves", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dx2 = params[1];
            num_t dy2 = params[2];
            num_t b = params[3];
            num_t e = params[4];
            num_t x = TX * b * sin(TY * dx2);
            num_t y = TY * e * sin(TX * dy2);
            VAR_RET(W * VEC(x,y));
        },
        // store: weight,dx2,dy2,b,e
        VAR_PARSE
        {
            const Affine2D<num_t>& aff = xform.getPreAffine();
            varp.push_back(weight);
            varp.push_back(1.0 / (aff.c*aff.c + EPS));
            varp.push_back(1.0 / (aff.f*aff.f + EPS));
            varp.push_back(aff.b);
            varp.push_back(aff.e);
        }
    )},
    {"fisheye", VAR_T(
        VAR_FUNC
        {
            num_t Wt2 = params[0]; // 2*weight
            num_t r = Wt2 / (TP.r() + 1.0);
            VAR_RET(r * VEC(TY,TX));
        },
        // store: 2*weight
        VAR_PARSE
        {
            varp.push_back(2.0*weight);
        }
    )},
    {"popcorn", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t c = params[1];
            num_t f = params[2];
            num_t dx = c * sin(tan(3.0*TY));
            num_t dy = f * sin(tan(3.0*TX));
            VAR_RET(W * (TP + VEC(dx,dy)));
        },
        // store: weight,c,f
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(xform.getPreAffine().c);
            varp.push_back(xform.getPreAffine().f);
        }
    )},
    {"exponential", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dx = W * exp(TX - 1.0);
            num_t sdy,cdy;
            sincosg(M_PI*TY,&sdy,&cdy);
            VAR_RET(dx * VEC(cdy,sdy));
        }
    )},
    {"power", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t sa,ca;
            sincosg(TP.atan(),&sa,&ca);
            num_t r = W * pow(TP.r(),sa);
            VAR_RET(r * VEC(ca,sa));
        }
    )},
    {"cosine", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t sa,ca;
            sincosg(TX*M_PI,&sa,&ca);
            VAR_RET(W * VEC(ca*cosh(TY),-sa*sinh(TY)));
        }
    )},
    {"rings", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dx = params[1];
            num_t r = TP.r();
            r = W * (fmod(r+dx,2.0*dx) - dx + r*(1.0-dx));
            num_t sa,ca;
            sincosg(TP.atan(),&sa,&ca);
            VAR_RET(r * VEC(ca,sa));
        },
        // store: weight,dx
        VAR_PARSE
        {
            const Affine2D<num_t>& aff = xform.getPreAffine();
            varp.push_back(weight);
            varp.push_back(aff.c*aff.c + eps<num_t>::value);
        }
    )},
    {"fan", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dx = params[1];
            num_t dy = params[2];
            static const num_t table[2] = {1.0,-1.0};
            num_t dx2 = dx*0.5;
            num_t a = TP.atan();
            num_t r = W * TP.r();
            num_t sa,ca;
            num_t m = table[fmod(a+dy,dx) > dx2];
            a += m*dx2;
            sincosg(a,&sa,&ca);
            VAR_RET(r * VEC(ca,sa));
        },
        // store: weight,dx,dy
        VAR_PARSE
        {
            const Affine2D<num_t>& aff = xform.getPreAffine();
            varp.push_back(weight);
            varp.push_back(M_PI*(aff.c*aff.c + EPS));
            varp.push_back(aff.f);
        }
    )},
    {"blob", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t mid = params[1];
            num_t diff2 = params[2];
            num_t waves = params[3];
            num_t r = TP.r();
            num_t a = TP.atan();
            r *= (mid + diff2*sin(waves*a));
            VAR_RET(W * r * VEC(sin(a),cos(a)));
        },
        // parse: low,high,waves
        // store: weight,mid,diff/2,waves
        VAR_PARSE
        {
            varp.push_back(weight);
            num_t low = json["low"].floatValue();
            num_t high = json["high"].floatValue();
            num_t waves = json["waves"].floatValue();
            varp.push_back((low+high)/2.0);
            varp.push_back((high-low)/2.0);
            varp.push_back(waves);
        }
    )},
    {"pdj", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = params[1];
            num_t b = params[2];
            num_t c = params[3];
            num_t d = params[4];
            num_t nx1 = cos(b*TX);
            num_t nx2 = sin(c*TX);
            num_t ny1 = sin(a*TY);
            num_t ny2 = cos(d*TY);
            VAR_RET(W * VEC(ny1-nx1,nx2-ny2));
        },
        // parse: a,b,c,d
        // store: weight,a,b,c,d
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["a"].floatValue());
            varp.push_back(json["b"].floatValue());
            varp.push_back(json["c"].floatValue());
            varp.push_back(json["d"].floatValue());
        }
    )},
    {"fan2", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dx = params[1];
            num_t dy = params[2];
            num_t dxinv = params[3];
            static const num_t table[2] = {1.0,-1.0};
            num_t dx2 = 0.5 * dx;
            num_t a = TP.atan();
            num_t sa,ca;
            num_t r = W * TP.r();
            num_t t = a + dy - dx*(i32)((a + dy) * dxinv);
            a += table[t > dx2]*dx2;
            sincosg(a,&sa,&ca);
            VAR_RET(r * VEC(sa,ca));
        },
        // parse: x,y
        // store: weight,dx,dy,1/dx
        VAR_PARSE
        {
            num_t x = json["x"].floatValue();
            num_t y = json["y"].floatValue();
            num_t dx = M_PI * (x*x + EPS);
            varp.push_back(weight);
            varp.push_back(dx);
            varp.push_back(y);
            varp.push_back(1.0/dx);
        }
    )},
    {"rings2", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dx = params[1];
            num_t dx2inv = params[2];
            num_t r = TP.r();
            r += -2.0*dx*(i32)((r+dx)*dx2inv) + r*(1.0 - dx);
            num_t a = TP.atan();
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(W * r * VEC(sa,ca));
        },
        // parse: value
        // store: weight,dx,1/(2*dx)
        VAR_PARSE
        {
            num_t v = json["value"].floatValue();
            num_t dx = v*v + EPS;
            varp.push_back(weight);
            varp.push_back(dx);
            varp.push_back(0.5/dx);
        }
    )},
    {"eyefish", VAR_T(
        VAR_FUNC
        {
            num_t Wt2 = params[0]; // 2*weight
            num_t r = Wt2 / (TP.r() + 1.0);
            VAR_RET(r * TP);
        },
        // store: 2*weight
        VAR_PARSE
        {
            varp.push_back(2.0*weight);
        }
    )},
    {"bubble", VAR_T(
        VAR_FUNC
        {
            num_t Wt4 = params[0];
            num_t r = Wt4 / (TP.r2() + 4.0); // W/(0.25*r^2+1) = 4*W/(r^2+4)
            VAR_RET(r * TP);
        },
        // store: 4*weight
        VAR_PARSE
        {
            varp.push_back(4.0*weight);
        }
    )},
    {"cylinder", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            VAR_RET(W * VEC(sin(TX),TY));
        }
    )},
    {"perspective", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t dist = params[1];
            num_t vsin = params[2];
            num_t vfcos = params[3];
            num_t t = 1.0 / (dist - TY*vsin);
            VAR_RET(W * t * VEC(dist*TX,vfcos*TY));
        },
        // parse: distance,angle
        // store: weight,distance,persp_vsin,persp_vfcos
        VAR_PARSE
        {
            num_t dist = json["distance"].floatValue();
            num_t angle = json["angle"].floatValue();
            varp.push_back(weight);
            varp.push_back(dist);
            varp.push_back(sin(angle));
            varp.push_back(dist*cos(angle));
        }
    )},
    {"noise", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t tr = (2.0 * M_PI) * state.randNum();
            num_t sr,cr;
            sincosg(tr,&sr,&cr);
            num_t r = W * state.randNum();
            VAR_RET(r * VEC(TX*cr,TY*sr));
        }
    )},
    {"julian", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t abspower = params[1];
            num_t invpower = params[2];
            num_t cn = params[3];
            i32 t = trunc(abspower*state.randNum());
            num_t a = (TP.atanyx() + (2.0*M_PI)*t) * invpower;
            num_t r = W * pow(TP.r2(),cn);
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(r * VEC(ca,sa));
        },
        // parse: power,distance
        // store: weight,abs(power),1/power,cn
        VAR_PARSE
        {
            num_t power = json["power"].floatValue();
            num_t dist = json["distance"].floatValue();
            varp.push_back(weight);
            varp.push_back(fabs(power));
            varp.push_back(1.0/power);
            varp.push_back(dist/(2.0*power));
        }
    )},
    {"juliascope", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t abspower = params[1];
            num_t invpower = params[2];
            num_t cn = params[3];
            static const num_t table[2] = {-1.0,1.0};
            i32 t = trunc(abspower*state.randNum());
            num_t a = ((2.0*M_PI)*t + table[t&1]*TP.atanyx()) * invpower;
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            num_t r = W * pow(TP.r2(),cn);
            VAR_RET(r * VEC(ca,sa));
        },
        // parse: power,distance
        // store: weight,abs(power),1/power,cn
        VAR_PARSE
        {
            num_t power = json["power"].floatValue();
            num_t dist = json["distance"].floatValue();
            varp.push_back(weight);
            varp.push_back(fabs(power));
            varp.push_back(1.0/power);
            varp.push_back(dist/(2.0*power));
        }
    )},
    {"blur", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = (2.0 * M_PI) * state.randNum();
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(W * state.randNum() * VEC(ca,sa));
        }
    )},
    {"gaussian_blur", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = (2.0 * M_PI) * state.randNum();
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            num_t g = (state.randNum() + state.randNum()
                + state.randNum() + state.randNum() - 2.0);
            VAR_RET(W * g * VEC(ca,sa));
        }
    )},
    {"radial_blur", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t spin = params[1];
            num_t zoom = params[2];
            num_t g = W * (state.randNum() + state.randNum()
                + state.randNum() + state.randNum() - 2.0);
            num_t ra = TP.r();
            num_t a = TP.atanyx() + spin*g;
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            num_t rz = zoom*g - 1.0;
            VAR_RET(ra*VEC(ca,sa) + rz*TP);
        },
        // parse: angle
        // store: weight,spin,zoom
        VAR_PARSE
        {
            num_t angle = json["angle"].floatValue();
            num_t spin,zoom;
            sincosg(angle*M_PI_2,&spin,&zoom);
            varp.push_back(weight);
            varp.push_back(spin);
            varp.push_back(zoom);
        }
    )},
    {"pie", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t slices = params[1];
            num_t rot = params[2];
            num_t thick = params[3];
            num_t invslices = params[4];
            i32 sl = (i32)(state.randNum()*slices + 0.5);
            num_t a = rot + (2.0*M_PI)*(sl + state.randNum()*thick)*invslices;
            num_t r = W * state.randNum();
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(r * VEC(ca,sa));
        },
        // parse: slices,rotation,thickness
        // store: weight,slices,rotation,thickness,1/slices
        VAR_PARSE
        {
            num_t slices = json["slices"].floatValue();
            varp.push_back(weight);
            varp.push_back(slices);
            varp.push_back(json["rotation"].floatValue());
            varp.push_back(json["thickness"].floatValue());
            varp.push_back(1.0/slices);
        }
    )},
    {"ngon", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t powerval = params[1];
            num_t angle = params[2];
            num_t corners = params[3];
            num_t circle = params[4];
            num_t invangle = params[5];
            static const num_t table[2] = {0.0,1.0};
            num_t r = pow(TP.r2(),powerval);
            num_t theta = TP.atanyx();
            num_t phi = theta - angle*floor(theta*invangle);
            phi -= table[phi > angle*0.5]*angle;
            num_t amp = corners*(1.0/(cos(phi)+EPS) - 1.0) + circle;
            amp /= (r + EPS);
            VAR_RET(W * amp * TP);
        },
        // parse: power,sides,corners,circle
        // store: weight,power/2,2*pi/sides,corners,circle,sides/(2*pi)
        VAR_PARSE
        {
            num_t sides = json["sides"].floatValue();
            varp.push_back(weight);
            varp.push_back(json["power"].floatValue()/2.0);
            varp.push_back(2.0*M_PI/sides);
            varp.push_back(json["corners"].floatValue());
            varp.push_back(json["circle"].floatValue());
            varp.push_back(sides/(2.0*M_PI));
        }
    )},
    {"curl", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t c1 = params[1];
            num_t c2 = params[2];
            num_t re = 1.0 + c1*TX + c2*(TX*TX - TY*TY);
            num_t im = c1*TY + 2.0*c2*TX*TY;
            num_t r = W / (re*re + im*im); // +EPS ???
            VAR_RET(r * VEC(TX*re+TY*im,TY*re-TX*im));
        },
        // parse: c1,c2
        // store: weight,c1,c2
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["c1"].floatValue());
            varp.push_back(json["c2"].floatValue());
        }
    )},
    {"rectangles", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t px = params[1];
            num_t py = params[2];
            num_t x,y;
            // branching with if/else seems to be the only good way here
            if (px == 0.0) x = TX;
            else x = (2.0*floor(TX/px) + 1.0)*px - TX;
            if (py == 0.0) y = TY;
            else y = (2.0*floor(TY/py) + 1.0)*py - TY;
            VAR_RET(W * VEC(x,y));
        },
        // parse: x,y
        // store: weight,x,y
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["x"].floatValue());
            varp.push_back(json["y"].floatValue());
        }
    )},
    {"arch", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = state.randNum() * W * M_PI;
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(W * VEC(sa,sa*sa/ca));
        }
    )},
    {"tangent", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            VAR_RET(W * VEC(sin(TX)/cos(TY),tan(TY)));
        }
    )},
    {"square", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t x = state.randNum();
            num_t y = state.randNum();
            VAR_RET(W * VEC(x-0.5,y-0.5));
        }
    )},
    {"rays", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t a = W * state.randNum() * M_PI;
            num_t r = W / (TP.r2() + EPS);
            num_t tr = W * tan(a) * r;
            VAR_RET(tr * VEC(cos(TX),sin(TY)));
        }
    )},
    {"blade", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r = state.randNum() * W * TP.r();
            num_t sr,cr;
            sincosg(r,&sr,&cr);
            VAR_RET(W * TX * VEC(cr+sr,cr-sr));
        }
    )},
    {"secant2", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            static const num_t table[2] = {-1.0,1.0};
            num_t cr = cos(W*TP.r());
            num_t icr = 1.0/cos(W*TP.r());
            VAR_RET(W * VEC(TX,icr+table[cr<0]));
        }
    )},
    {"twintrian", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r = state.randNum() * W * TP.r();
            num_t sr,cr,diff;
            sincosg(r,&sr,&cr);
            diff = log10(sr*sr)+cr;
            if (unlikely(bad_value(diff))) diff = -30.0;
            VAR_RET(W * TX * VEC(diff,diff-sr*M_PI));
        }
    )},
    {"cross", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s = TX*TX - TY*TY;
            num_t r = W * sqrt(1.0 / (s*s + EPS));
            VAR_RET(r * TP);
        }
    )},
    {"disc2", VAR_T(
        VAR_FUNC
        {
            num_t W_pi = params[0];
            num_t rotpi = params[1];
            num_t cosadd = params[2];
            num_t sinadd = params[3];
            num_t t = rotpi * (TX + TY);
            num_t st,ct;
            sincosg(t,&st,&ct);
            num_t r = W_pi * TP.atan();
            VAR_RET(r * VEC(st+cosadd,ct+sinadd));
        },
        // parse: rotation,twist
        // store: weight/pi,timespi,cosadd,sinadd
        VAR_PARSE
        {
            varp.push_back(weight*M_1_PI);
            num_t rot = json["rotation"].floatValue();
            num_t twist = json["twist"].floatValue();
            varp.push_back(rot*M_PI);
            num_t cosadd,sinadd;
            sincosg(twist,&sinadd,&cosadd);
            cosadd -= 1.0;
            num_t k = 1.0;
            if (twist > 2.0*M_PI)
                k = (1.0+twist-2.0*M_PI);
            if (twist < -2.0*M_PI)
                k = (1.0+twist+2.0*M_PI);
            varp.push_back(cosadd*k);
            varp.push_back(sinadd*k);
        }
    )},
    {"supershape", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t pm_4 = params[1];
            num_t pneg1_n1 = params[2];
            num_t n2 = params[3];
            num_t n3 = params[4];
            num_t rnd = params[5];
            num_t holes = params[6];
            num_t theta = pm_4*TP.atanyx() + M_PI_4;
            num_t st,ct;
            sincosg(theta,&st,&ct);
            num_t t1 = pow(fabs(ct),n2);
            num_t t2 = pow(fabs(st),n3);
            num_t tr = TP.r();
            num_t r = W * ((rnd*state.randNum() + (1.0-rnd)*tr) - holes)
                        * pow(t1+t2,pneg1_n1) / tr; // +EPS ???
            VAR_RET(r * TP);
        },
        // parse: rnd,m,n1,n2,n3,holes
        // store: weight,pm_4,pneg1_n1,n2,n3,rnd,holes
        VAR_PARSE
        {
            varp.push_back(weight);
            num_t n1 = json["n1"].floatValue();
            varp.push_back(json["m"].floatValue()/4.0);
            varp.push_back(-1.0/n1);
            varp.push_back(json["n2"].floatValue());
            varp.push_back(json["n3"].floatValue());
            varp.push_back(json["rnd"].floatValue());
            varp.push_back(json["holes"].floatValue());
        }
    )},
    {"flower", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t petals = params[1];
            num_t holes = params[2];
            num_t theta = TP.atanyx();
            num_t r = W * (state.randNum() - holes) * cos(petals*theta)
                / TP.r(); // +EPS ???
            VAR_RET(r * TP);
        },
        // parse: petals,holes
        // store: weight,petals,holes
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["petals"].floatValue());
            varp.push_back(json["holes"].floatValue());
        }
    )},
    {"conic", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t eccen = params[1];
            num_t holes = params[2];
            num_t tr = TP.r(); // +EPS ???
            num_t ct = TX / tr;
            num_t r = W * (state.randNum() - holes) * eccen
                / (tr + tr*eccen*ct);
            VAR_RET(r * TP);
        },
        // parse: eccentricity,holes
        // store: weight,eccentricity,holes
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["eccentricity"].floatValue());
            varp.push_back(json["holes"].floatValue());
        }
    )},
    {"parabola", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t ht = params[1];
            num_t wt = params[2];
            num_t sr,cr;
            sincosg(TP.r(),&sr,&cr);
            VAR_RET(W * VEC(ht*sr*sr*state.randNum(),wt*cr*state.randNum()));
        },
        // parse: height,width
        // store: weight,height,width
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["height"].floatValue());
            varp.push_back(json["width"].floatValue());
        }
    )},
    {"bent2", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t px = params[1];
            num_t py = params[2];
            // TODO eliminate if statements
            num_t nx = TX;
            num_t ny = TY;
            if (nx < 0.0) nx *= px;
            if (ny < 0.0) ny *= py;
            VAR_RET(W * VEC(nx,ny));
        },
        // parse: x,y
        // store: weight,x,y
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["x"].floatValue());
            varp.push_back(json["y"].floatValue());
        }
    )},
    {"bipolar", VAR_T(
        VAR_FUNC
        {
            num_t W2_pi = params[0];
            num_t shift = params[1];
            num_t x2y2 = TP.r2();
            num_t t = x2y2+1.0;
            num_t x2 = 2.0*TX;
            num_t y = 0.5*atan2(2.0*TY,x2y2-1.0) + shift;
            y -= M_PI * floor(y*M_1_PI + 0.5);
            VAR_RET(W2_pi * VEC(0.25*log((t+x2)/(t-x2)),y));
        },
        // parse: shift
        // store: weight*2/pi,-(pi/2)*shift
        VAR_PARSE
        {
            varp.push_back(weight*M_2_PI);
            varp.push_back(-M_PI_2*json["shift"].floatValue());
        }
    )},
    {"boarders", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            static const num_t table[2] = {-0.25,0.25};
            // TODO ensure rounding works right, seems default is nearest
            num_t rx = rint(TX);
            num_t ry = rint(TY);
            num_t ox = TX - rx;
            num_t oy = TY - ry;
            if (state.randNum() >= 0.75)
                VAR_RET(W * VEC(ox*0.5+rx,oy*0.5+ry));
            else
            {
                // TODO eliminate some branches
                if (fabs(ox) >= fabs(oy))
                {
                    bool z = ox >= 0.0;
                    num_t x = ox*0.5+rx+table[z];
                    num_t y = oy*0.5+ry+table[z]*oy/ox;
                    VAR_RET(W * VEC(x,y));
                }
                else
                {
                    bool z = oy >= 0.0;
                    num_t x = ox*0.5+rx+table[z]*ox/oy;
                    num_t y = oy*0.5+ry+table[z];
                    VAR_RET(W * VEC(x,y));
                }
            }
        }
    )},
    {"butterfly", VAR_T(
        VAR_FUNC
        {
            num_t wx = params[0];
            num_t y2 = 2.0*TY;
            num_t r = wx * sqrt(fabs(TX*TY)  / (EPS + TX*TX + y2*y2));
            VAR_RET(r * VEC(TX,y2));
        },
        // store: weight*1.3029400317411197908970256609023
        VAR_PARSE
        {
            // constant 4/sqrt(3*pi) from flam3 source
            varp.push_back(weight*1.3029400317411197908970256609023);
        }
    )},
    {"cell", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t size = params[1];
            num_t invsize = params[2];
            num_t x = floor(TX * invsize);
            num_t y = floor(TY * invsize);
            num_t dx = TX - x*size;
            num_t dy = TY - y*size;
            // TODO can branches be eliminated
            if (y >= 0.0)
            {
                if (x >= 0.0) x *= 2.0, y *= 2.0;
                else y *= 2.0, x = -(2.0*x+1.0);
            }
            else
            {
                if (x >= 0.0) y = -(2.0*y+1.0), x *= 2.0;
                else y = -(2.0*y+1.0), x = -(2.0*x+1.0);
            }
            VAR_RET(W * VEC(dx+x*size,-dy-y*size));
        },
        // parse: size
        // store: weight,size,1/size
        VAR_PARSE
        {
            num_t size = json["size"].floatValue();
            varp.push_back(weight);
            varp.push_back(size);
            varp.push_back(1.0/size);
        }
    )},
    {"cpow", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t va = params[1];
            num_t vc = params[2];
            num_t vd = params[3];
            num_t power = params[4];
            num_t a = TP.atanyx();
            num_t lnr = 0.5 * log(TP.r2());
            num_t ang = vc*a + vd*lnr + va*floor(power*state.randNum());
            num_t sa,ca;
            sincosg(ang,&sa,&ca);
            num_t m = W * exp(vc*lnr - vd*a);
            VAR_RET(m * VEC(ca,sa));
        },
        // parse: r,i,power
        // store: weight,va,vc,vd,power
        VAR_PARSE
        {
            num_t r = json["r"].floatValue();
            num_t i = json["i"].floatValue();
            num_t power = json["power"].floatValue();
            varp.push_back(weight);
            varp.push_back(2.0*M_PI/power);
            varp.push_back(r/power);
            varp.push_back(i/power);
            varp.push_back(power);
        }
    )},
    {"curve", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t invxl = params[1];
            num_t invyl = params[2];
            num_t xamp = params[3];
            num_t yamp = params[4];
            VEC_T v = VEC(xamp*exp(-TY*TY*invxl),yamp*(-TX*TX*invyl));
            VAR_RET(W * (TP + v));
        },
        // parse: xamp,yamp,xlen,ylen
        // store: weight,1.0/pc_xlen,1.0/pc_ylen,xamp,yamp
        VAR_PARSE
        {
            num_t xamp = json["xamp"].floatValue();
            num_t yamp = json["yamp"].floatValue();
            num_t xlen = json["xlen"].floatValue();
            num_t ylen = json["ylen"].floatValue();
            num_t pc_xlen = xlen*xlen;
            num_t pc_ylen = ylen*ylen;
            varp.push_back(weight);
            varp.push_back(1.0 / (pc_xlen < 1e-20 ? 1e-20 : pc_xlen));
            varp.push_back(1.0 / (pc_ylen < 1e-20 ? 1e-20 : pc_ylen));
            varp.push_back(xamp);
            varp.push_back(yamp);
        }
    )},
    {"edisc", VAR_T(
        VAR_FUNC
        {
            num_t W_c = params[0];
            static const num_t table[2] = {1.0,-1.0};
            num_t tmp = TP.r2() + 1.0;
            num_t tmp2 = 2.0 * TX;
            num_t xmax = 0.5*(sqrt(tmp+tmp2)+sqrt(tmp-tmp2));
            num_t a1 = log(xmax + sqrt(xmax - 1.0));
            num_t a2 = -acos(TX / xmax);
            num_t s1,c1;
            sincosg(a1,&s1,&c1);
            num_t s2 = sinh(a2);
            num_t c2 = cosh(a2);
            s1 *= table[TY > 0.0];
            VAR_RET(W_c * VEC(c2*c1,s2*s1));
        },
        // store: weight*0.0864278365005759 (division by 11.57034632 in flam3)
        VAR_PARSE
        {
            varp.push_back(weight*0.0864278365005759);
        }
    )},
    {"elliptic", VAR_T(
        VAR_FUNC
        {
            num_t W2_pi = params[0];
            static const num_t table[2] = {-1.0,1.0};
            num_t tmp = TP.r2() + 1.0;
            num_t x2 = 2.0 * TX;
            num_t xmax = 0.5 * (sqrt(tmp+x2)+sqrt(tmp-x2));
            num_t a = TX / xmax;
            num_t b = 1.0 - a*a;
            num_t ssx = xmax - 1.0;
            // TODO can branches be eliminated
            b = b < 0.0 ? 0.0 : sqrt(b);
            ssx = ssx < 0.0 ? 0.0 : sqrt(ssx);
            VAR_RET(W2_pi * VEC(atan2(a,b),table[TY>0.0]*log(xmax+ssx)));
        },
        // store: weight*2/pi
        VAR_PARSE
        {
            varp.push_back(weight*M_2_PI);
        }
    )},
    {"escher", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t vc = params[1];
            num_t vd = params[2];
            num_t a = TP.atanyx();
            num_t lnr = 0.5 * log(TP.r2());
            num_t m = W * exp(vc*lnr - vd*a);
            num_t n = vc*a + vd*lnr;
            num_t sn,cn;
            sincosg(n,&sn,&cn);
            VAR_RET(m * VEC(cn,sn));
        },
        // parse: beta
        // store: weight,vc,vd
        VAR_PARSE
        {
            num_t beta = json["beta"].floatValue();
            num_t seb,ceb;
            sincosg(beta,&seb,&ceb);
            varp.push_back(weight);
            varp.push_back(0.5*(1.0+ceb));
            varp.push_back(0.5*seb);
        }
    )},
    {"foci", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t expx = 0.5 * exp(TX);
            num_t expnx = 0.25 / expx;
            num_t sn,cn;
            sincosg(TY,&sn,&cn);
            num_t tmp = W / (expx + expnx - cn);
            VAR_RET(tmp * VEC(expx-expnx,sn));
        }
    )},
    {"lazysusan", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t px = params[1];
            num_t py = params[2];
            num_t spin = params[3];
            num_t twist = params[4];
            num_t spacep1tw = params[5];
            num_t x = TX - px;
            num_t y = TY - py;
            num_t r = hypot(x,y); // +EPS ???
            num_t sa,ca;
            if (r < W)
            {
                num_t a = atan2(y,x) + spin + twist*(W-r);
                sincosg(a,&sa,&ca);
                r *= W;
                VAR_RET(VEC(r*ca+px,r*sa-py));
            }
            else
            {
                r = spacep1tw / r;
                VAR_RET(VEC(r*x+px,r*y-py));
            }
        },
        // parse: spin,space,twist,x,y
        // store: weight,x,y,spin,twist,(1+space)*weight
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["x"].floatValue());
            varp.push_back(json["y"].floatValue());
            varp.push_back(json["spin"].floatValue());
            varp.push_back(json["twist"].floatValue());
            varp.push_back((1.0+json["space"].floatValue())*weight);
        }
    )},
    {"loonie", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t r2 = TP.r2(); // +EPS ???
            num_t w2 = W*W;
            num_t r = W;
            if (r2 < w2) r *= sqrt(w2/r2 - 1.0);
            VAR_RET(r * TP);
        }
    )},
    {"pre_blur", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t g = W * (state.randNum() + state.randNum()
                + state.randNum() + state.randNum() - 2.0);
            num_t a = (2.0 * M_PI) * state.randNum();
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(g * VEC(ca,sa));
        }
    )},
    {"modulus", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t px = params[1];
            num_t py = params[2];
            num_t pxinv = params[3];
            num_t pyinv = params[4];
            num_t x = TX - px*floor(TX*pxinv + 0.5);
            num_t y = TY - py*floor(TY*pyinv + 0.5);
            VAR_RET(W * VEC(x,y));
        },
        // parse: x,y
        // store: weight,2*x,2*y,1/(2*x),1/(2*y)
        VAR_PARSE
        {
            num_t x = json["x"].floatValue();
            num_t y = json["y"].floatValue();
            varp.push_back(weight);
            varp.push_back(2.0*x);
            varp.push_back(2.0*y);
            varp.push_back(1.0/(2.0*x));
            varp.push_back(1.0/(2.0*y));
        }
    )},
    {"oscope", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t tpf = params[1];
            num_t p_amp = params[2];
            num_t p_damp = params[3];
            num_t sep = params[4];
            static const num_t table[2] = {1.0,-1.0};
            num_t damp = exp(-fabs(TX)*p_damp);
            num_t t = p_amp * damp * cos(tpf*TX) + sep;
            num_t y = table[fabs(TY) <= t] * TY;
            VAR_RET(W * VEC(TX,y));
        },
        // parse: separation,frequency,amplitude,damping
        // store: weight,tpf(2*pi*freq),amplitude,damping,separation
        VAR_PARSE
        {
            num_t freq = json["frequency"].floatValue();
            varp.push_back(weight);
            varp.push_back(2.0*M_PI*freq);
            varp.push_back(json["amplitude"].floatValue());
            varp.push_back(json["damping"].floatValue());
            varp.push_back(json["separation"].floatValue());
        }
    )},
    {"polar2", VAR_T(
        VAR_FUNC
        {
            num_t W_pi = params[0];
            VAR_RET(W_pi * VEC(TP.atan(),0.5*log(TP.r2())));
        },
        // store: weight/pi
        VAR_PARSE
        {
            varp.push_back(weight*M_1_PI);
        }
    )},
    {"popcorn2", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t px = params[1];
            num_t py = params[2];
            num_t pc = params[3];
            num_t dx = px*sin(tan(TY*pc));
            num_t dy = py*sin(tan(TX*pc));
            VAR_RET(W * (TP + VEC(dx,dy)));
        },
        // parse: x,y,c
        // store: weight,x,y,c
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["x"].floatValue());
            varp.push_back(json["y"].floatValue());
            varp.push_back(json["c"].floatValue());
        }
    )},
    {"scry", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t t = TP.r2();
            num_t r = 1.0 / (sqrt(t) * (t + 1.0/(W + EPS)));
            VAR_RET(r * TP);
        }
    )},
    {"separation", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t x2 = params[1];
            num_t y2 = params[2];
            num_t xin = params[3];
            num_t yin = params[4];
            static const num_t table[2] = {-1.0,1.0};
            bool xp = TX > 0.0;
            bool yp = TY > 0.0;
            num_t x = sqrt(TX*TX + x2) + table[!xp]*xin;
            num_t y = sqrt(TY*TY + y2) + table[!yp]*yin;
            VAR_RET(W * VEC(table[xp]*x,table[yp]*y));
        },
        // parse: x,y,xin,yin
        // store: weight,x*x,y*y,xin,yin
        VAR_PARSE
        {
            num_t x = json["x"].floatValue();
            num_t y = json["y"].floatValue();
            varp.push_back(weight);
            varp.push_back(x*x);
            varp.push_back(y*y);
            varp.push_back(json["xin"].floatValue());
            varp.push_back(json["yin"].floatValue());
        }
    )},
    {"split", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t xspi = params[1];
            num_t yspi = params[2];
            static const num_t table[2] = {-1.0,1.0};
            bool xp = cos(TX*xspi) >= 0.0;
            bool yp = cos(TY*yspi) >= 0.0;
            VAR_RET(W * VEC(table[xp]*TY,table[yp]*TX));
        },
        // parse: xsize,ysize
        // store: weight,pi*xsize,pi*ysize
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(M_PI*json["xsize"].floatValue());
            varp.push_back(M_PI*json["ysize"].floatValue());
        }
    )},
    {"splits", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t px = params[1];
            num_t py = params[2];
            static const num_t table[2] = {-1.0,1.0};
            VAR_RET(W * (TP + VEC(table[TX>=0.0]*px,table[TY>=0.0]*py)));
        },
        // parse: x,y
        // store: weight,x,y
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["x"].floatValue());
            varp.push_back(json["y"].floatValue());
        }
    )},
    {"stripes", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t space = params[1];
            num_t warp = params[2];
            num_t rx = floor(TX + 0.5);
            num_t ox = TX - rx;
            VAR_RET(W * VEC(ox*space+rx,TY+ox*ox*warp));
        },
        // parse: space,warp
        // store: weight,1-space,warp
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(1.0-json["space"].floatValue());
            varp.push_back(json["warp"].floatValue());
        }
    )},
    {"wedge", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t swirl = params[1];
            num_t count = params[2];
            num_t angle = params[3];
            num_t hole = params[4];
            num_t r = TP.r();
            num_t a = TP.atanyx() + swirl*r;
            num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
            num_t cf = 1.0 - angle*swirl*(M_1_PI*0.5);
            a = a*cf + c*angle;
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(W * (r + hole) * VEC(ca,sa));
        },
        // parse: angle,hole,count,swirl
        // store: weight,swirl,count,angle,hole
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["swirl"].floatValue());
            varp.push_back(json["count"].floatValue());
            varp.push_back(json["angle"].floatValue());
            varp.push_back(json["hole"].floatValue());
        }
    )},
    {"wedge_julia", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t cn = params[1];
            num_t abspower = params[2];
            num_t power = params[3];
            num_t count = params[4];
            num_t angle = params[5];
            num_t cf = params[6];
            num_t r = W * pow(TP.r2(),cn);
            i32 tr = (i32)(abspower * state.randNum());
            num_t a = (TP.atanyx() + (2.0*M_PI)*tr) / power;
            num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
            num_t sa,ca;
            a = a*cf + c*angle;
            sincosg(a,&sa,&ca);
            VAR_RET(r * VEC(ca,sa));
        },
        // parse: angle,count,power,distance
        // store: weight,cn,abs(power),power,count,angle,cf
        VAR_PARSE
        {
            num_t angle = json["angle"].floatValue();
            num_t count = json["count"].floatValue();
            num_t power = json["power"].floatValue();
            num_t dist = json["distance"].floatValue();
            varp.push_back(weight);
            varp.push_back(dist/(2.0*power));
            varp.push_back(fabs(power));
            varp.push_back(power);
            varp.push_back(count);
            varp.push_back(angle);
            varp.push_back(1.0-angle*count*M_1_PI*0.5);
        }
    )},
    {"wedge_sph", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t swirl = params[1];
            num_t count = params[2];
            num_t cf = params[3];
            num_t angle = params[4];
            num_t hole = params[5];
            num_t r = 1.0 / (TP.r() + EPS);
            num_t a = TP.atanyx() + swirl*r;
            num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
            num_t sa,ca;
            a = a*cf + c*angle;
            sincosg(a,&sa,&ca);
            VAR_RET(W * (r + hole) * VEC(ca,sa));
        },
        // parse: angle,count,hole,swirl
        // store: weight,swirl,count,cf,angle,hole
        VAR_PARSE
        {
            num_t angle = json["angle"].floatValue();
            num_t count = json["count"].floatValue();
            varp.push_back(weight);
            varp.push_back(json["swirl"].floatValue());
            varp.push_back(count);
            varp.push_back(1.0-angle*count*M_1_PI*0.5);
            varp.push_back(angle);
            varp.push_back(json["hole"].floatValue());
        }
    )},
    {"whorl", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t choice[2] = {params[1],params[2]};
            num_t r = TP.r();
            num_t a = TP.atanyx();
            a += choice[r >= W] / (W - r);
            num_t sa,ca;
            sincosg(a,&sa,&ca);
            VAR_RET(W * r * VEC(ca,sa));
        },
        // parse: inside,outside
        // store: weight,inside,outside
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["inside"].floatValue());
            varp.push_back(json["outside"].floatValue());
        }
    )},
    {"waves2", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t xfreq = params[1];
            num_t xscale = params[2];
            num_t yfreq = params[3];
            num_t yscale = params[4];
            num_t dx = xscale*sin(TY*xfreq);
            num_t dy = yscale*sin(TX*yfreq);
            VAR_RET(W * (TP + VEC(dx,dy)));
        },
        // parse: xfreq,xscale,yfreq,yscale
        // store: weight,xfreq,xscale,yfreq,yscale
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["xfreq"].floatValue());
            varp.push_back(json["xscale"].floatValue());
            varp.push_back(json["yfreq"].floatValue());
            varp.push_back(json["yscale"].floatValue());
        }
    )},
    {"exp", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t e = exp(TX);
            num_t es,ec;
            sincosg(TY,&es,&ec);
            VAR_RET(W * e * VEC(ec,es));
        }
    )},
    {"log", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            VAR_RET(W * VEC(0.5*log(TP.r2()),TP.atanyx()));
        }
    )},
    {"sin", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(TX,&s,&c);
            num_t sh = sinh(TY);
            num_t ch = cosh(TY);
            VAR_RET(W * VEC(s*ch,c*sh));
        }
    )},
    {"cos", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(TX,&s,&c);
            num_t sh = sinh(TY);
            num_t ch = cosh(TY);
            VAR_RET(W * VEC(c*ch,-s*sh));
        }
    )},
    {"tan", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(2.0*TX,&s,&c);
            num_t sh = sinh(2.0*TY);
            num_t ch = cosh(2.0*TY);
            num_t den = W/(c+ch);
            VAR_RET(den * VEC(s,sh));
        }
    )},
    {"sec", VAR_T(
        VAR_FUNC
        {
            num_t Wt2 = params[0];
            num_t s,c;
            sincosg(TX,&s,&c);
            num_t sh = sinh(TY);
            num_t ch = cosh(TY);
            num_t den = Wt2/(cos(2.0*TX)+cosh(2.0*TY));
            VAR_RET(den * VEC(c*ch,s*sh));
        },
        // store: 2*weight
        VAR_PARSE
        {
            varp.push_back(2.0*weight);
        }
    )},
    {"csc", VAR_T(
        VAR_FUNC
        {
            num_t Wt2 = params[0];
            num_t s,c;
            sincosg(TX,&s,&c);
            num_t sh = sinh(TY);
            num_t ch = cosh(TY);
            num_t den = Wt2/(cosh(2.0*TY)-cos(2.0*TX));
            VAR_RET(den * VEC(s*ch,-c*sh));
        },
        // store: 2*weight
        VAR_PARSE
        {
            varp.push_back(2.0*weight);
        }
    )},
    {"cot", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(2.0*TX,&s,&c);
            num_t sh = sinh(2.0*TY);
            num_t ch = cosh(2.0*TY);
            num_t den = W/(ch-c);
            VAR_RET(den * VEC(s,-sh));
        }
    )},/////////////////////////////////////////////////////////////////////////
    {"sinh", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(TY,&s,&c);
            num_t sh = sinh(TX);
            num_t ch = cosh(TX);
            VAR_RET(W * VEC(sh*c,ch*s));
        }
    )},
    {"cosh", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(TY,&s,&c);
            num_t sh = sinh(TX);
            num_t ch = cosh(TX);
            VAR_RET(W * VEC(ch*c,sh*s));
        }
    )},
    {"tanh", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(2.0*TY,&s,&c);
            num_t sh = sinh(2.0*TX);
            num_t ch = cosh(2.0*TX);
            num_t den = W/(c+ch);
            VAR_RET(den * VEC(sh,s));
        }
    )},
    {"sech", VAR_T(
        VAR_FUNC
        {
            num_t Wt2 = params[0];
            num_t s,c;
            sincosg(TY,&s,&c);
            num_t sh = sinh(TX);
            num_t ch = cosh(TX);
            num_t den = Wt2/(cos(2.0*TY)+cosh(2.0*TX));
            VAR_RET(den * VEC(c*ch,-s*sh));
        },
        // store: 2*weight
        VAR_PARSE
        {
            varp.push_back(2.0*weight);
        }
    )},
    {"csch", VAR_T(
        VAR_FUNC
        {
            num_t Wt2 = params[0];
            num_t s,c;
            sincosg(TY,&s,&c);
            num_t sh = sinh(TX);
            num_t ch = cosh(TX);
            num_t den = Wt2/(cosh(2.0*TX)-cos(2.0*TY));
            VAR_RET(den * VEC(sh*c,-ch*s));
        },
        // store: 2*weight
        VAR_PARSE
        {
            varp.push_back(2.0*weight);
        }
    )},
    {"coth", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t s,c;
            sincosg(2.0*TY,&s,&c);
            num_t sh = sinh(2.0*TX);
            num_t ch = cosh(2.0*TX);
            num_t den = W/(ch-c);
            VAR_RET(den * VEC(sh,s));
        }
    )},
    {"auger", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t freq = params[1];
            num_t aw = params[2];
            num_t scale = params[3];
            num_t sym = params[4];
            num_t s = sin(freq*TX);
            num_t t = sin(freq*TY);
            num_t dy = TY + aw*(scale + fabs(TY))*s;
            num_t dx = aw*(scale + fabs(TX))*t;
            VAR_RET(W * VEC(TX+sym*dx,dy));
        },
        // parse: sym,auger_weight,freq,scale
        // store: weight,freq,auger_weight,scale/2,sym
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["freq"].floatValue());
            varp.push_back(json["auger_weight"].floatValue());
            varp.push_back(json["scale"].floatValue()/2.0);
            varp.push_back(json["sym"].floatValue());
        }
    )},
    {"flux", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t spread = params[1];
            num_t xpw = TX + W;
            num_t xmw = TX - W;
            num_t y2 = TY * TY;
            num_t ar = W * spread * sqrt(sqrt(y2 + xpw*xpw)/sqrt(y2 + xmw*xmw));
            num_t aa = (atan2(TY,xmw) - atan2(TY,xpw)) * 0.5;
            num_t sa,ca;
            sincosg(aa,&sa,&ca);
            VAR_RET(ar * VEC(ca,sa));
        },
        // parse: spread
        // store: weight,2+spread
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(2.0+json["spread"].floatValue());
        }
    )},
    {"mobius", VAR_T(
        VAR_FUNC
        {
            num_t W = params[0];
            num_t re_a = params[1];
            num_t re_b = params[2];
            num_t re_c = params[3];
            num_t re_d = params[4];
            num_t im_a = params[5];
            num_t im_b = params[6];
            num_t im_c = params[7];
            num_t im_d = params[8];
            num_t re_u = re_a*TX - im_a*TY + re_b;
            num_t im_u = re_a*TY + im_a*TX + im_b;
            num_t re_v = re_c*TX - im_c*TY + re_d;
            num_t im_v = re_c*TY + im_c*TX + im_d;
            num_t rad_v = W / (re_v*re_v + im_v*im_v); // +EPS ???
            VAR_RET(rad_v * VEC(re_u*re_v+im_u*im_v,im_u*re_v-re_u*im_v));
        },
        // parse: re_a,re_b,re_c,re_d,im_a,im_b,im_c,im_d
        // store: weight,re_a,re_b,re_c,re_d,im_a,im_b,im_c,im_d
        VAR_PARSE
        {
            varp.push_back(weight);
            varp.push_back(json["re_a"].floatValue());
            varp.push_back(json["re_b"].floatValue());
            varp.push_back(json["re_c"].floatValue());
            varp.push_back(json["re_d"].floatValue());
            varp.push_back(json["im_a"].floatValue());
            varp.push_back(json["im_b"].floatValue());
            varp.push_back(json["im_c"].floatValue());
            varp.push_back(json["im_d"].floatValue());
        }
    )}
};

}
}

#undef likely
#undef unlikely

#undef VAR_T
#undef STATE_T
#undef XFORM_T
#undef JSON_T
#undef PARAM_T
#undef VAR_FUNC
#undef VAR_PARSE
#undef VAR_RET
#undef TX
#undef TY
#undef TP
#undef EPS
