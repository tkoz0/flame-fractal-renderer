/*
Variation functions

TODO separate into multiple files
TODO better error messages
*/

#pragma once

#include "../rng/flame_rng.hpp"

#include "../types/point.hpp"
#include "../types/types.hpp"
#include "../types/constants.hpp"
#include "../types/flame.hpp"

#include "../utils/json.hpp"
#include "../utils/math.hpp"
#include "../utils/sfinae.hpp"
#include "../utils/flame.hpp"

#include <cstdlib>
#include <string>

// macros for using SFINAE
#define ENABLE_IF(COND,RET) template <typename RET2 = RET> \
    [[nodiscard]] std::enable_if_t<(COND),RET2>
#define ENABLE_IF2(COND) template <typename RET = Variation<dims>*> \
    [[nodiscard]] std::enable_if_t<(COND),RET>
#define ENABLE_IF3(COND) template <typename RET> \
    [[nodiscard]] std::enable_if_t<(COND),RET>

// forward declarations
namespace tkoz::flame::vars
{

// variation base class
template <size_t dims>
class Variation
{
private:
    num_t weight;
    Variation(){} // must provide json for constructor
public:
    Variation(const Json& json)
    {
        weight = json["weight"].floatValue();
    }
    virtual ~Variation(){}
    [[nodiscard]] virtual inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const = 0;
    [[nodiscard]] inline num_t getWeight() const { return weight; }
    ENABLE_IF2(dims<2) static parseVariation2d(const Json& json);
    ENABLE_IF2(dims>=2) static parseVariation2d(const Json& json);
    ENABLE_IF2(dims<3) static parseVariation3d(const Json& json);
    ENABLE_IF2(dims>=3) static parseVariation3d(const Json& json);
    [[nodiscard]] static Variation<dims> *parseVariationNd(const Json& json);
    [[nodiscard]] static Variation<dims> *parseVariation(const Json& json);
};

// class for generalizing 2d variation to higher dimensions
// specify axes to operate on with index keys axis_x and axis_y
template <size_t dims>
class VariationFrom2D: public Variation<dims>
{
    static_assert(dims >= 2);
private:
    u32 axis_x,axis_y;
    typedef Point<num_t,dims> _nd;
    typedef Point<num_t,2> _2d;
public:
    VariationFrom2D(const Json& json): Variation<dims>(json)
    {
        if (dims < 2)
            throw std::runtime_error("minimum of 2 dimensions required");
        if (dims == 2)
            return;
        JsonInt ax = json["axis_x"].intValue();
        JsonInt ay = json["axis_y"].intValue();
        if (ax < 0 || ax >= (JsonInt)dims)
            throw std::runtime_error("axis_x index out of range");
        if (ay < 0 || ay >= (JsonInt)dims)
            throw std::runtime_error("axis_y index out of range");
        if (ax == ay)
            throw std::runtime_error("axes are not distinct");
        axis_x = ax;
        axis_y = ay;
    }
    virtual ~VariationFrom2D(){}
    [[nodiscard]] inline _nd calc(rng_t& rng, const _nd& tx) const
    {
        return calc_h(rng,tx);
    }
    ENABLE_IF(dims==2,_nd) inline calc_h(
        rng_t& rng, const _nd& tx) const
    {
        return calc2d(rng,tx);
    }
    ENABLE_IF(dims>2,_nd) inline calc_h(
        rng_t& rng, const _nd& tx) const
    {
        _2d ret2d = calc2d(rng,_2d(tx[axis_x],tx[axis_y]));
        _nd ret;
        ret[axis_x] = ret2d.x();
        ret[axis_y] = ret2d.y();
        return ret;
    }
    [[nodiscard]] virtual inline _2d calc2d(
        rng_t& rng, const _2d& tx) const = 0;
};

// class for generalizing 3d variation to higher dimensions
// specify axes to operate on with index keys axis_x and axis_y and axis_z
template <size_t dims>
class VariationFrom3D: public Variation<dims>
{
    static_assert(dims >= 3);
private:
    u32 axis_x,axis_y,axis_z;
    typedef Point<num_t,dims> _nd;
    typedef Point<num_t,3> _3d;
public:
    VariationFrom3D(const Json& json): Variation<dims>(json)
    {
        if (dims < 3)
            throw std::runtime_error("minimum of 2 dimensions required");
        if (dims == 3)
            return;
        JsonInt ax = json["axis_x"].intValue();
        JsonInt ay = json["axis_y"].intValue();
        JsonInt az = json["axis_z"].intValue();
        if (ax < 0 || ax >= (JsonInt)dims)
            throw std::runtime_error("axis_x index out of range");
        if (ay < 0 || ay >= (JsonInt)dims)
            throw std::runtime_error("axis_y index out of range");
        if (az < 0 || az >= (JsonInt)dims)
            throw std::runtime_error("axis_z index out of range");
        if (ax == ay || ax == az || ay == az)
            throw std::runtime_error("axes are not distinct");
        axis_x = ax;
        axis_y = ay;
        axis_z = az;
    }
    virtual ~VariationFrom3D(){}
    [[nodiscard]] inline _nd calc(
        rng_t& rng, const _nd& tx) const
    {
        return calc_h(rng,tx);
    }
    ENABLE_IF(dims==3,_nd) inline calc_h(
        rng_t& rng, const _nd& tx) const
    {
        return calc2d(rng,tx);
    }
    ENABLE_IF(dims>3,_nd) inline calc_h(
        rng_t& rng, const _nd& tx) const
    {
        _3d ret3d = calc3d(rng,_3d(tx[axis_x],tx[axis_y],tx[axis_z]));
        _nd ret;
        ret[axis_x] = ret3d.x();
        ret[axis_y] = ret3d.y();
        ret[axis_z] = ret3d.z();
        return ret;
    }
    [[nodiscard]] virtual inline _3d calc3d(
        rng_t& rng, const _3d& tx) const = 0;
};

/*
======= N dimensional variations from flam3 =======
*/

/*
Linear - generalized from flam3 var0_linear
*/
template <size_t dims>
struct Linear: public Variation<dims>
{
    Linear(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        return tx;
    }
};

/*
Sinusoidal - generalized from flam3 var1_sinusoidal
Replace each coordinate with its sine
*/
template <size_t dims>
struct Sinusoidal: public Variation<dims>
{
    Sinusoidal(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        return tx.map([](num_t x){ return sin(x); });
    }
};

/*
Spherical - generalized from flam3 var2_spherical
Divide by squared 2-norm
*/
template <size_t dims>
struct Spherical: public Variation<dims>
{
    Spherical(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        // for small dimensions maybe division
        // is better than precomputing 1/r^2
        num_t r = 1.0 / (tx.norm2sq() + eps_v<num_t>);
        return r * tx;
    }
};

/*
Bent - generalized from flam3 var14_bent and var54_bent2
In flam3, the default would be neg[2.0,0.5],pos[1,1] for bent
The bent2 variation allows scaling negatives but not positives
*/
template <size_t dims>
class Bent: public Variation<dims>
{
    Point<num_t,dims> scales_neg;
    Point<num_t,dims> scales_pos;
public:
    Bent(const Json& json): Variation<dims>(json)
    {
        // should these default to all ones
        scales_neg = Point<num_t,dims>(json["scales_neg"]);
        scales_pos = Point<num_t,dims>(json["scales_pos"]);
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret = tx;
        for (size_t i = 0; i < dims; ++i)
        {
            if (ret[i] < 0.0)
                ret[i] *= scales_neg[i];
            else // flam3 does not scale positive coordinates
                ret[i] *= scales_pos[i];
        }
        return ret;
    }
};

/*
Rectangles - generalized from flam3 var40_rectangles
*/
template <size_t dims>
class Rectangles: public Variation<dims>
{
    Point<num_t,dims> params;
public:
    Rectangles(const Json& json): Variation<dims>(json)
    {
        params = Point<num_t,dims>(json["params"]);
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret;
        for (size_t i = 0; i < dims; ++i)
        {
            num_t p = params[i];
            num_t x = tx[i];
            if (p == 0.0)
                ret[i] = x;
            else
                ret[i] = (2.0*floor(x/p) + 1.0)*p - x;
        }
        return ret;
    }
};

/*
Fisheye - generalized from flam3 var16_fisheye and var27_eyefish
Keeping correct order from eyefish, rather than incorrect order of fisheye
The original behavior can be achieved with affine transformations
*/
template <size_t dims>
class Fisheye: public Variation<dims>
{
    num_t addval;
public:
    Fisheye(const Json& json): Variation<dims>(json)
    {
        // should this be required to be positive?
        // in flam3, this is 1, should this be default?
        addval = json["addval"].floatValue();
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        // is division faster than precomputing 1/(r+addval)
        num_t r = 1.0 / (tx.norm2() + addval);
        return r * tx;
    }
};

/*
Bubble - generalized from flam3 var28_bubble
*/
template <size_t dims>
class Bubble: public Variation<dims>
{
    num_t addval;
public:
    Bubble(const Json& json): Variation<dims>(json)
    {
        // should this be required positive or default to 4?
        // default is 4 in flam3
        addval = json["addval"].floatValue();
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        // is division faster or is precomputing 1/(r^2+addval) faster?
        num_t r = 1.0 / (tx.norm2sq() + addval);
        return r * tx;
    }
};

/*
Noise - generalized from flam3 var31_noise
*/
template <size_t dims>
struct Noise: public Variation<dims>
{
    Noise(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        num_t r = rng.randNum();
        // flam3 uses random unit vector in 2d plane
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return r * multComponents(tx,dir);
    }
};

/*
Blur - generalized from flam3 var34_blur
*/
template <size_t dims>
struct Blur: public Variation<dims>
{
    Blur(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        num_t r = rng.randNum();
        // flam3 uses random unit vector in 2d plane
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return r * dir;
    }
};

/*
Gaussian Blur - generalized from flam3 var35_gaussian
*/
template <size_t dims>
struct GaussianBlur: public Variation<dims>
{
    GaussianBlur(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        // flam3 simulates gaussian distribution with (4 randoms [0,1)) - 2
        num_t r = rng.randGaussian();
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return r * dir;
    }
};

/*
Square Noise - generalized from flam3 var43_square
*/
template <size_t dims>
struct SquareNoise: public Variation<dims>
{
    SquareNoise(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        return rng.template randPoint2<dims>();
    }
};

/*
Separation - generalized from flam3 var73_separation
*/
template <size_t dims>
class Separation: public Variation<dims>
{
    Point<num_t,dims> sep2,inside;
public:
    Separation(const Json& json): Variation<dims>(json)
    {
        sep2 = Point<num_t,dims>(json["params"]);
        for (size_t i = 0; i < dims; ++i)
            sep2[i] *= sep2[i];
        inside = Point<num_t,dims>(json["inside"]);
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret;
        for (size_t i = 0; i < dims; ++i)
        {
            num_t s = copysign(1.0,tx[i]);
            ret[i] = s * (sqrt(tx[i]*tx[i] + sep2[i]) - s*inside[i]);
        }
        return ret;
    }
};

/*
Splits - generalized from flam3 var75_splits
*/
template <size_t dims>
class Splits: public Variation<dims>
{
    Point<num_t,dims> params;
public:
    Splits(const Json& json): Variation<dims>(json)
    {
        params = Point<num_t,dims>(json["params"]);
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret;
        for (size_t i = 0; i < dims; ++i)
        {
            num_t s = copysign(1.0,tx[i]);
            ret[i] = s*params[i];
        }
        return tx + ret;
    }
};

/*
Pre Blur - generalized from flam3 var67_pre_blur
*/
template <size_t dims>
struct PreBlur: public Variation<dims>
{
    PreBlur(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        // flam3 simulates gaussian with (sum 4 randoms in [0,1)) - 2
        num_t g = rng.randGaussian();
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return g * dir;
    }
};

/*
Modulus - generalized from flam3 var68_modulus
*/
template <size_t dims>
class Modulus: public Variation<dims>
{
    Point<num_t,dims> params2, params2inv;
public:
    Modulus(const Json& json): Variation<dims>(json)
    {
        params2 = Point<num_t,dims>(json["params"]);
        params2 *= 2.0;
        params2inv = params2.map([](num_t x){ return 1.0/x; });
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret;
        // computed in a different way, should behave the same as flam3
        for (size_t i = 0; i < dims; ++i)
            ret[i] = tx[i] - params2[i]*floor(tx[i]*params2inv[i] + 0.5);
        return ret;
    }
};

/*
Cell N - generalized from flam3 var58_cell
Flam3 uses the same size for both dimensions
*/
template <size_t dims>
class CellN: public Variation<dims>
{
    Point<num_t,dims> sizes,invsizes;
public:
    CellN(const Json& json): Variation<dims>(json)
    {
        sizes = Point<num_t,dims>(json["sizes"]);
        for (size_t i = 0; i < dims; ++i)
            invsizes[i] = 1.0 / sizes[i];
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret;
        for (size_t i = 0; i < dims; ++i)
        {
            num_t x = floor(tx[i] * invsizes[i]);
            num_t dx = tx[i] - x*sizes[i];
            num_t xs = copysign(2.0,x);
            num_t x2 = x * xs;
            x2 -= (num_t)(x < 0);
            ret[i] = dx + x2*sizes[i];
        }
        return ret;
    }
};

/*
======= 2 dimensional variations from flam3 =======
*/

/*
Swirl - 2d, from flam3 var3_swirl
*/
template <size_t dims>
struct Swirl: public VariationFrom2D<dims>
{
    Swirl(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = tx.norm2sq();
        num_t x,y;
        tx.getXY(x,y);
        num_t sr,cr;
        math::sincosg(r,sr,cr);
        return Point<num_t,2>(x*sr-y*cr,x*cr+y*sr);
    }
};

/*
Horseshoe - 2d, from flam3 var4_horseshoe
*/
template <size_t dims>
struct Horseshoe: public VariationFrom2D<dims>
{
    Horseshoe(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        // is division better than precomputing 1/r
        num_t r = 1.0 / (tx.norm2() + eps_v<num_t>);
        num_t x,y;
        tx.getXY(x,y);
        // is it faster to compute x^2-y^2
        return r * Point<num_t,2>((x-y)*(x+y),2.0*x*y);
    }
};

/*
Polar - 2d, from flam3 var5_polar
*/
template <size_t dims>
struct Polar: public VariationFrom2D<dims>
{
    Polar(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        return Point<num_t,2>(a*M_1_PI,r-1.0);
    }
};

/*
Polar2 - 2d, from flam3 var70_polar2
*/
template <size_t dims>
struct Polar2: public VariationFrom2D<dims>
{
    Polar2(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        // leaving out 0.5 factor on the log
        return Point<num_t,2>(tx.angle(),log(tx.norm2sq()));
    }
};

/*
Handkerchief - 2d, from flam3 var6_handkerchief
*/
template <size_t dims>
struct Handkerchief: public VariationFrom2D<dims>
{
    Handkerchief(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        return r * Point<num_t,2>(sin(a+r),cos(a-r));
    }
};

/*
Heart - 2d, from flam3 var7_heart
*/
template <size_t dims>
struct Heart: public VariationFrom2D<dims>
{
    Heart(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        num_t sa,ca;
        math::sincosg(r*a,sa,ca);
        return r * Point<num_t,2>(sa,-ca);
    }
};

/*
Disc - 2d, from flam3 var8_disc
*/
template <size_t dims>
struct Disc: public VariationFrom2D<dims>
{
    Disc(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle(); // not including M_1_PI factor
        num_t r = tx.norm2();
        num_t sr,cr;
        math::sincosg(M_PI*r,sr,cr);
        return a * Point<num_t,2>(sr,cr);
    }
};

/*
Disc2 - 2d, from flam3 var49_disc2
*/
template <size_t dims>
class Disc2: public VariationFrom2D<dims>
{
    num_t rotpi;
    Point<num_t,2> addval;
public:
    Disc2(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t rot = json["rotation"].floatValue();
        num_t twist = json["twist"].floatValue();
        rotpi = rot*M_PI;
        num_t sinadd,cosadd;
        math::sincosg(twist,sinadd,cosadd);
        cosadd -= 1.0;
        num_t k = 1.0 + twist;
        if (twist > 2.0*M_PI) k -= 2.0*M_PI;
        if (twist < -2.0*M_PI) k += 2.0*M_PI;
        cosadd *= k;
        sinadd *= k;
        addval = Point<num_t,2>(cosadd,sinadd);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t t = rotpi * (x + y);
        num_t st,ct;
        math::sincosg(t,st,ct);
        // order of ct,st switched from flam3
        return tx.angle() * (Point<num_t,2>(ct,st) + addval);
    }
};

/*
Waves - 2d, from flam3 var15_waves and var81_waves2
*/
template <size_t dims>
class Waves: public VariationFrom2D<dims>
{
    num_t xf,xs,yf,ys;
public:
    Waves(const Json& json): VariationFrom2D<dims>(json)
    {
        xf = json["xfreq"].floatValue();
        xs = json["xscale"].floatValue();
        yf = json["yfreq"].floatValue();
        ys = json["yscale"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t dx = xs*sin(y*xf);
        num_t dy = ys*sin(x*yf);
        return tx + Point<num_t,2>(dx,dy);
    }
};

/*
Fan - 2d, from flam3 var22_fan and var25_fan2
*/
template <size_t dims>
class Fan: public VariationFrom2D<dims>
{
    num_t dx,dy;
public:
    Fan(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t x = json["x"].floatValue();
        num_t y = json["y"].floatValue();
        dx = M_PI * (x*x + eps_v<num_t>);
        dy = y;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t dx2 = dx*0.5;
        num_t a = tx.angle();
        // the fmod part should behave the same as t in flam3
        num_t m = copysign(1.0,dx2-fmod(a+dy,dx));
        a += m*dx2;
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        // order of ca,sa switched from flam3
        return tx.norm2() * Point<num_t,2>(ca,sa);
    }
};

/*
Rings - 2d, from flam3 var21_rings and var26_rings2
*/
template <size_t dims>
class Rings: public VariationFrom2D<dims>
{
    num_t dx;
public:
    Rings(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t v = json["value"].floatValue();
        dx = v*v + eps_v<num_t>;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,s,c;
        tx.getRadiusSinCos(r,s,c);
        // should behave the same as in flam3
        r = fmod(r+dx,2.0*dx) - dx + r*(1.0-dx);
        // order of c,s switched from flam3
        return r * Point<num_t,2>(c,s);
    }
};

/*
Spiral - 2d, from flam3 var9_spiral
*/
template <size_t dims>
struct Spiral: public VariationFrom2D<dims>
{
    Spiral(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        num_t sr,cr;
        math::sincosg(r,sr,cr);
        // is division faster than precomputing 1/r
        num_t r1 = 1.0 / (r + eps_v<num_t>);
        return r1 * Point<num_t,2>(ca+sr,sa-cr);
    }
};

/*
Hyperbolic - 2d, from flam3 var10_hyperbolic
*/
template <size_t dims>
struct Hyperbolic: public VariationFrom2D<dims>
{
    Hyperbolic(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        return Point<num_t,2>(sa/(r+eps_v<num_t>),ca*r);
    }
};

/*
Diamond - 2d, from flam3 var11_diamond
*/
template <size_t dims>
struct Diamond: public VariationFrom2D<dims>
{
    Diamond(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        num_t sr,cr;
        math::sincosg(r,sr,cr);
        return Point<num_t,2>(sa*cr,ca*sr);
    }
};

/*
Ex - 2d, from flam3 var12_ex
*/
template <size_t dims>
struct Ex: public VariationFrom2D<dims>
{
    Ex(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        num_t n0 = sin(a+r);
        num_t n1 = cos(a-r);
        num_t m0 = n0*n0*n0 * r;
        num_t m1 = n1*n1*n1 * r;
        return Point<num_t,2>(m0+m1,m0-m1);
    }
};

/*
Julia - 2d, from flam3 var13_julia
*/
template <size_t dims>
struct Julia: public VariationFrom2D<dims>
{
    Julia(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t a = 0.5*tx.angle() + rng.randBool()*M_PI; // + 0 or PI
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        // switched order of sin,cos from flam3
        return tx.norm2() * Point<num_t,2>(ca,sa);
    }
};

/*
Exponential - 2d, from flam3 var18_exponential
*/
template <size_t dims>
struct Exponential: public VariationFrom2D<dims>
{
    Exponential(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t dx = exp(x-1.0);
        num_t sdy,cdy;
        math::sincosg(M_PI*y,sdy,cdy);
        return dx * Point<num_t,2>(cdy,sdy);
    }
};

/*
Power - 2d, from flam3 var19_power
*/
template <size_t dims>
struct Power: public VariationFrom2D<dims>
{
    Power(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        return pow(r,sa) * Point<num_t,2>(ca,sa);
    }
};

/*
Cosine - 2d, from flam3 var20_cosine
*/
template <size_t dims>
struct Cosine: public VariationFrom2D<dims>
{
    Cosine(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t sa,ca;
        math::sincosg(x*M_PI,sa,ca);
        return Point<num_t,2>(ca*cosh(y),-sa*sinh(y));
    }
};

/*
Blob - 2d, from flam3 var23_blob
*/
template <size_t dims>
class Blob: public VariationFrom2D<dims>
{
    num_t mid,amp,waves;
public:
    Blob(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t low = json["low"].floatValue();
        num_t high = json["high"].floatValue();
        mid = (high+low)/2.0;
        amp = (high-low)/2.0;
        waves = json["waves"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        num_t a = tx.angle();
        r *= mid + amp*sin(waves*a); // simplified
        return r * Point<num_t,2>(ca,sa); // order switched
    }
};

/*
PDJ - 2d, from flam3 var24_pdj
*/
template <size_t dims>
class PDJ: public VariationFrom2D<dims>
{
    num_t a,b,c,d;
public:
    PDJ(const Json& json): VariationFrom2D<dims>(json)
    {
        a = json["a"].floatValue();
        b = json["b"].floatValue();
        c = json["c"].floatValue();
        d = json["d"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t nx1 = cos(b*x);
        num_t nx2 = sin(c*x);
        num_t ny1 = sin(a*y);
        num_t ny2 = cos(d*y);
        return Point<num_t,2>(ny1-nx1,nx2-ny2);
    }
};

/*
Cylinder - 2d, from flam3 var29_cylinder
Should this variation be replaced with something more general?
Similar to sinusoidal but only applies on some axes
*/
template <size_t dims>
struct Cylinder: public VariationFrom2D<dims>
{
    Cylinder(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        return Point<num_t,2>(sin(tx.x()),tx.y());
    }
};

/*
Perspective - 2d, from flam3 var30_perspective
*/
template <size_t dims>
class Perspective: public VariationFrom2D<dims>
{
    num_t dist,vsin,vfcos;
public:
    Perspective(const Json& json): VariationFrom2D<dims>(json)
    {
        dist = json["distance"].floatValue();
        num_t angle = json["angle"].floatValue();
        vsin = sin(angle);
        vfcos = dist*cos(angle);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t t = 1.0 / (dist - y*vsin);
        return t * Point<num_t,2>(dist*x,vfcos*y);
    }
};

/*
Julia N - 2d, from flam3 var32_juliaN_generic
*/
template <size_t dims>
class JuliaN: public VariationFrom2D<dims>
{
    num_t abspower,invpower,cn;
public:
    JuliaN(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t power = json["power"].floatValue();
        num_t dist = json["dist"].floatValue();
        abspower = fabs(power); // flam3 calls this rN
        invpower = 1.0/power;
        cn = dist/(2.0*power);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        i32 t = trunc(abspower*rng.randNum());
        num_t a = (tx.angle() + (2.0*M_PI)*t) * invpower;
        num_t r = pow(tx.norm2sq(),cn);
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Julia Scope - 2d, from flam3 var33_juliaScope_generic
*/
template <size_t dims>
class JuliaScope: public VariationFrom2D<dims>
{
    num_t abspower,invpower,cn;
public:
    JuliaScope(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t power = json["power"].floatValue();
        num_t dist = json["dist"].floatValue();
        abspower = fabs(power); // flam3 calls this rN
        invpower = 1.0/power;
        cn = dist/(2.0*power);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        i32 t = trunc(abspower*rng.randNum());
        num_t dir = rng.template randDirection<1>().x();
        num_t a = ((2.0*M_PI)*t + dir*tx.angle()) * invpower;
        num_t r = pow(tx.norm2sq(),cn);
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Radial Blur - 2d, from flam3 var36_radial_blur
Flam3 uses weight in a nonstandard way
*/
template <size_t dims>
class RadialBlur: public VariationFrom2D<dims>
{
    num_t spin,zoom,flam3weight;
public:
    RadialBlur(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t angle = json["angle"].floatValue();
        math::sincosg(angle*M_PI_2,spin,zoom);
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        // flam3 simulates gaussian with (4 randoms in [0,1)) - 2
        num_t g = flam3weight * rng.randGaussian();
        num_t ra = tx.norm2();
        num_t a = tx.angle() + spin*g;
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        num_t rz = zoom*g - 1.0;
        return ra*Point<num_t,2>(ca,sa) + rz*tx;
    }
};

/*
Pie - 2d, from flam3 var37_pie
*/
template <size_t dims>
class Pie: public VariationFrom2D<dims>
{
    num_t slices,rotation,thickness,invslices2pi;
public:
    Pie(const Json& json): VariationFrom2D<dims>(json)
    {
        slices = json["slices"].floatValue();
        rotation = json["rotation"].floatValue();
        thickness = json["thickness"].floatValue();
        invslices2pi = (2.0*M_PI)/slices;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)tx;
        i32 sl = (i32)(rng.randNum()*slices + 0.5);
        num_t a = rotation + (sl + rng.randNum()*thickness)*invslices2pi;
        num_t r = rng.randNum();
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
NGon - 2d, from flam3 var38_ngon
*/
template <size_t dims>
class NGon: public VariationFrom2D<dims>
{
    num_t powerval,angle,corners,circle,invangle;
public:
    NGon(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t sides = json["sides"].floatValue();
        powerval = json["power"].floatValue()/2.0;
        angle = (2.0*M_PI)/sides; // in flam3, named b, not precomputed
        corners = json["corners"].floatValue();
        circle = json["circle"].floatValue();
        invangle = sides/(2.0*M_PI);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = pow(tx.norm2sq(),powerval);
        num_t theta = tx.angle();
        num_t phi = theta - angle*floor(theta*invangle);
        static const num_t mult[2] = {0.0,1.0};
        phi -= mult[phi > angle*0.5]*angle;
        num_t amp = corners*(1.0/(cos(phi)+eps_v<num_t>)-1.0) + circle;
        amp /= r + eps_v<num_t>;
        return amp * tx;
    }
};

/*
Curl - 2d, from flam3 var39_curl
*/
template <size_t dims>
class Curl: public VariationFrom2D<dims>
{
    num_t c1,c2;
public:
    Curl(const Json& json): VariationFrom2D<dims>(json)
    {
        c1 = json["c1"].floatValue();
        c2 = json["c2"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t re = 1.0 + c1*x + c2*(x*x - y*y);
        num_t im = c1*y + 2.0*c2*x*y;
        // is it faster to divide than precompute 1/(re^2+im^2)
        // flam3 does not add eps to the denominator
        num_t r = 1.0 / (re*re + im*im + eps_v<num_t>);
        return r * Point<num_t,2>(x*re+y*im,y*re-x*im);
    }
};

/*
Arch - 2d, from flam3 var41_arch
Uses a nonstandard weight
*/
template <size_t dims>
class Arch: public VariationFrom2D<dims>
{
    num_t flam3weight;
public:
    Arch(const Json& json): VariationFrom2D<dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)tx;
        num_t a = flam3weight * rng.randNum() * M_PI;
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        return Point<num_t,2>(sa,sa*sa/ca);
    }
};

/*
Tangent - 2d, from flam3 var42_tangent
*/
template <size_t dims>
struct Tangent: public VariationFrom2D<dims>
{
    Tangent(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        return Point<num_t,2>(sin(x)/cos(y),tan(y));
    }
};

/*
Rays - 2d, from flam3 var44_rays
Uses nonstandard variation weight
*/
template <size_t dims>
class Rays: public VariationFrom2D<dims>
{
    num_t flam3weight;
public:
    Rays(const Json& json): VariationFrom2D<dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t a = flam3weight * rng.randNum() * M_PI;
        num_t r = flam3weight / (tx.norm2sq() + eps_v<num_t>);
        num_t tr = tan(a) * r; // not multiplying flam3weight here
        return tr * Point<num_t,2>(cos(x),sin(y));
    }
};

/*
Blade - 2d, from flam3 var45_blade
Uses nonstandard variation weight
*/
template <size_t dims>
class Blade: public VariationFrom2D<dims>
{
    num_t flam3weight;
public:
    Blade(const Json& json): VariationFrom2D<dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t r = rng.randNum() * flam3weight * tx.norm2();
        num_t sr,cr;
        math::sincosg(r,sr,cr);
        return tx.x() * Point<num_t,2>(cr+sr,cr-sr);
    }
};

/*
Secant - 2d, from flam3 var46_secant2
Uses nonstandard variation weight
*/
template <size_t dims>
class Secant: public VariationFrom2D<dims>
{
    num_t flam3weight;
public:
    Secant(const Json& json): VariationFrom2D<dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t cr = cos(flam3weight*tx.norm2());
        num_t icr = 1.0/cr;
        num_t sign = copysign(1.0,-cr);
        return Point<num_t,2>(tx.x(),icr+sign);
    }
};

/*
Twintrian - 2d, from flam3 var47_twintrian
Uses nonstandard variation weight
*/
template <size_t dims>
class Twintrian: public VariationFrom2D<dims>
{
    num_t flam3weight;
public:
    Twintrian(const Json& json): VariationFrom2D<dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t r = rng.randNum() * flam3weight * tx.norm2();
        num_t sr,cr;
        math::sincosg(r,sr,cr);
        num_t diff = log10(sr*sr) + cr; // why log10?
        // flam3 uses this, does it make sense?
        if (bad_value(diff)) [[unlikely]]
            diff = -30.0;
        return tx.x() * Point<num_t,2>(diff,diff-sr*M_PI);
    }
};

/*
Cross - 2d, from flam3 var48_cross
*/
template <size_t dims>
struct Cross: public VariationFrom2D<dims>
{
    Cross(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s = x*x - y*y;
        num_t r = sqrt(1.0 / (s*s + eps_v<num_t>));
        return r * tx;
    }
};

/*
Exp - 2d, from flam3 var82_exp
*/
template <size_t dims>
struct Exp: public VariationFrom2D<dims>
{
    Exp(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t e = exp(x);
        num_t es,ec;
        math::sincosg(y,es,ec);
        return e * Point<num_t,2>(ec,es);
    }
};

/*
Log - 2d, from flam3 var83_log
*/
template <size_t dims>
struct Log: public VariationFrom2D<dims>
{
    Log(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        return Point<num_t,2>(log(tx.norm2sq()),tx.angle());
    }
};

/*
Sin - 2d, from flam3 var84_sin
*/
template <size_t dims>
struct Sin: public VariationFrom2D<dims>
{
    Sin(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(x,s,c);
        num_t sh = sinh(y);
        num_t ch = cosh(y);
        return Point<num_t,2>(s*ch,c*sh);
    }
};

/*
Cos - 2d, from flam3 var85_cos
*/
template <size_t dims>
struct Cos: public VariationFrom2D<dims>
{
    Cos(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(x,s,c);
        num_t ch = cosh(y);
        num_t sh = sinh(y);
        return Point<num_t,2>(c*ch,-s*sh);
    }
};

/*
Tan - 2d, from flam3 var86_tan
*/
template <size_t dims>
struct Tan: public VariationFrom2D<dims>
{
    Tan(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(2.0*x,s,c);
        num_t sh = sinh(2.0*y);
        num_t ch = cosh(2.0*y);
        // is it faster to precompute 1/(c+ch) and multiply?
        return Point<num_t,2>(s,sh) / (c+ch);
    }
};

/*
Sec - 2d, from flam3 var87_sec
*/
template <size_t dims>
struct Sec: public VariationFrom2D<dims>
{
    Sec(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(x,s,c);
        num_t sh = sinh(y);
        num_t ch = cosh(y);
        // is it faster to compute 1/(cos(2x)+cosh(2y)) and multiply?
        return Point<num_t,2>(c*ch,s*sh) / (cos(2.0*x)+cosh(2.0*y));
    }
};

/*
Csc - 2d, from flam3 var88_csc
*/
template <size_t dims>
struct Csc: public VariationFrom2D<dims>
{
    Csc(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(x,s,c);
        num_t sh = sinh(y);
        num_t ch = cosh(y);
        // is it faster to compute 1/(cosh(2y)-cos(2x)) and multiply?
        return Point<num_t,2>(s*ch,-c*sh) / (cosh(2.0*y)-cos(2.0*x));
    }
};

/*
Cot - 2d, from flam3 var89_cot
*/
template <size_t dims>
struct Cot: public VariationFrom2D<dims>
{
    Cot(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(2.0*x,s,c);
        num_t sh = sinh(2.0*y);
        num_t ch = cosh(2.0*y);
        // is it faster to compute 1/(ch-c) and multiply?
        return Point<num_t,2>(s,-sh) / (ch-c);
    }
};

/*
Sinh - 2d, from flam3 var90_sinh
*/
template <size_t dims>
struct Sinh: public VariationFrom2D<dims>
{
    Sinh(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        return Point<num_t,2>(sh*c,ch*s);
    }
};

/*
Cosh - 2d, from flam3 var91_cosh
*/
template <size_t dims>
struct Cosh: public VariationFrom2D<dims>
{
    Cosh(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        return Point<num_t,2>(ch*c,sh*s);
    }
};

/*
Tanh - 2d, from flam3 var92_tanh
*/
template <size_t dims>
struct Tanh: public VariationFrom2D<dims>
{
    Tanh(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(2.0*y,s,c);
        num_t sh = sinh(2.0*x);
        num_t ch = cosh(2.0*x);
        // is it faster to compute 1/(ch+c) and multiply?
        return Point<num_t,2>(sh,s) / (c+ch);
    }
};

/*
Sech - 2d, from flam3 var93_sech
*/
template <size_t dims>
struct Sech: public VariationFrom2D<dims>
{
    Sech(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        // is it faster to compute 1/(cos(2y)+cosh(2x)) and multiply?
        return Point<num_t,2>(c*ch,-s*sh) / (cos(2.0*y)+cosh(2.0*x));
    }
};

/*
Csch - 2d, from flam3 var94_csch
*/
template <size_t dims>
struct Csch: public VariationFrom2D<dims>
{
    Csch(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        // is it faster to compute 1/(cosh(2x)-cos(2y)) and multiply?
        return Point<num_t,2>(sh*c,-ch*s) / (cosh(2.0*x)-cos(2.0*y));
    }
};

/*
Coth - 2d, from flam3 var95_coth
*/
template <size_t dims>
struct Coth: public VariationFrom2D<dims>
{
    Coth(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        math::sincosg(2.0*y,s,c);
        num_t sh = sinh(2.0*x);
        num_t ch = cosh(2.0*x);
        // is it faster to compute 1/(ch-c) and multiply?
        return Point<num_t,2>(sh,s) / (ch-c);
    }
};

/*
Auger - 2d, from flam3 var96_auger
*/
template <size_t dims>
class Auger: public VariationFrom2D<dims>
{
    num_t freq,augerweight,scale,sym;
public:
    Auger(const Json& json): VariationFrom2D<dims>(json)
    {
        freq = json["freq"].floatValue();
        augerweight = json["flam3_weight"].floatValue();
        scale = json["scale"].floatValue() / 2.0;
        sym = json["sym"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s = sin(freq*x);
        num_t t = sin(freq*y);
        num_t dy = y + augerweight*(scale + fabs(y))*s;
        num_t dx = x + augerweight*(scale + fabs(x))*t;
        return Point<num_t,2>(x+sym*(dx-x),dy);
    }
};

/*
Flux - 2d, from flam3 var97_flux
*/
template <size_t dims>
class Flux: public VariationFrom2D<dims>
{
    num_t spread,fluxweight;
public:
    Flux(const Json& json): VariationFrom2D<dims>(json)
    {
        // why does flam3 have a +2 for this?
        spread = 2.0 + json["spread"].floatValue();
        fluxweight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t xpw = x + fluxweight;
        num_t xmw = x - fluxweight;
        num_t y2 = y*y;
        num_t avgr = spread * sqrt(sqrt(y2+xpw*xpw)/sqrt(y2+xmw*xmw));
        num_t avga = (atan2(y,xmw) - atan2(y,xpw)) * 0.5;
        num_t c,s;
        math::sincosg(avga,c,s);
        return avgr * Point<num_t,2>(c,s);
    }
};

/*
Mobius - 2d, from flam3 var98_mobius
*/
template <size_t dims>
class Mobius: public VariationFrom2D<dims>
{
    Point<num_t,2> a,b,c,d;
public:
    Mobius(const Json& json): VariationFrom2D<dims>(json)
    {
        a = Point<num_t,2>(json["a"]);
        b = Point<num_t,2>(json["b"]);
        c = Point<num_t,2>(json["c"]);
        d = Point<num_t,2>(json["d"]);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t re_u = a.x()*x - a.y()*y + b.x();
        num_t im_u = a.x()*y + a.y()*x + b.y();
        num_t re_v = c.x()*x - c.y()*y + d.x();
        num_t im_v = c.x()*y + c.y()*x + d.y();
        // faster to divide instead of precomputing this?
        num_t rad = 1.0 / (re_v*re_v + im_v*im_v + eps_v<num_t>);
        return rad * Point<num_t,2>(re_u*re_v+im_u*im_v,im_u*re_v-re_u*im_v);
    }
};

/*
Scry - 2d, from flam3 var72_scry
Uses a nonstandard variation weight
*/
template <size_t dims>
class Scry: public VariationFrom2D<dims>
{
    num_t scryweight;
public:
    Scry(const Json& json): VariationFrom2D<dims>(json)
    {
        scryweight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t t = tx.norm2sq();
        num_t r = 1.0 / (sqrt(t) * (t + 1.0/(scryweight + eps_v<num_t>)));
        return r * tx;
    }
};

/*
Split - 2d, from flam3 var74_split
*/
template <size_t dims>
class Split: public VariationFrom2D<dims>
{
    num_t sizex,sizey;
public:
    Split(const Json& json): VariationFrom2D<dims>(json)
    {
        sizex = json["xsize"].floatValue() * M_PI;
        sizey = json["ysize"].floatValue() * M_PI;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t xs = copysign(1.0,cos(x*sizex));
        num_t ys = copysign(1.0,cos(y*sizey));
        return Point<num_t,2>(x*ys,y*xs);
    }
};

/*
Stripes - 2d, from flam3 var76_stripes
*/
template <size_t dims>
class Stripes: public VariationFrom2D<dims>
{
    num_t space,warp;
public:
    Stripes(const Json& json): VariationFrom2D<dims>(json)
    {
        space = 1.0 - json["space"].floatValue();
        warp = json["warp"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t rx = floor(x + 0.5);
        num_t ox = x - rx;
        return Point<num_t,2>(ox*space+rx,y+ox*ox*warp);
    }
};

/*
Wedge - 2d, from flam3 var77_wedge
*/
template <size_t dims>
class Wedge: public VariationFrom2D<dims>
{
    num_t swirl,count,angle,hole,cf;
public:
    Wedge(const Json& json): VariationFrom2D<dims>(json)
    {
        swirl = json["swirl"].floatValue();
        count = json["count"].floatValue();
        angle = json["angle"].floatValue();
        hole = json["hole"].floatValue();
        cf = 1.0 - angle*count*(M_1_PI*0.5);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = tx.norm2();
        num_t a = tx.angle() + swirl*r;
        num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
        a = a*cf + c*angle;
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        return (r+hole) * Point<num_t,2>(ca,sa);
    }
};

/*
Wedge Julia - 2d, from flam3 var78_wedge_julia
*/
template <size_t dims>
class WedgeJulia: public VariationFrom2D<dims>
{
    num_t cn,abspower,invpower,count,angle,cf;
public:
    WedgeJulia(const Json& json): VariationFrom2D<dims>(json)
    {
        angle = json["angle"].floatValue();
        count = json["count"].floatValue();
        num_t power = json["power"].floatValue();
        invpower = 1.0/power;
        num_t dist = json["dist"].floatValue();
        cn = dist/(2.0*power);
        abspower = fabs(power); // flam3 calls this rN
        cf = 1.0 - angle*count*(M_1_PI*0.5);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t r = pow(tx.norm2sq(),cn);
        i32 tr = (i32)(abspower * rng.randNum());
        num_t a = (tx.angle() + (2.0*M_PI)*tr) * invpower;
        num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
        num_t sa,ca;
        a = a*cf + c*angle;
        math::sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Wedge Sph - 2d, from flam3 var79_wedge_sph
*/
template <size_t dims>
class WedgeSph: public VariationFrom2D<dims>
{
    num_t swirl,count,cf,angle,hole;
public:
    WedgeSph(const Json& json): VariationFrom2D<dims>(json)
    {
        angle = json["angle"].floatValue();
        count = json["count"].floatValue();
        swirl = json["swirl"].floatValue();
        hole = json["hole"].floatValue();
        cf = 1.0 - angle*count*(M_1_PI*0.5);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2() + eps_v<num_t>);
        num_t a = tx.angle() + swirl*r;
        num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
        num_t sa,ca;
        a = a*cf + c*angle;
        math::sincosg(a,sa,ca);
        return (r+hole) * Point<num_t,2>(ca,sa);
    }
};

/*
Whorl - 2d, from flam3 var80_whorl
Uses a nonstandard variation weight
*/
template <size_t dims>
class Whorl: public VariationFrom2D<dims>
{
    num_t choice[2],whorlweight;
public:
    Whorl(const Json& json): VariationFrom2D<dims>(json)
    {
        choice[0] = json["inside"].floatValue();
        choice[1] = json["outside"].floatValue();
        whorlweight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = tx.norm2();
        num_t a = tx.angle();
        a += choice[r >= whorlweight] / (whorlweight - r);
        num_t sa,ca;
        math::sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Supershape - 2d, from flam3 var50_supershape
*/
template <size_t dims>
class Supershape: public VariationFrom2D<dims>
{
    num_t pm_4,pneg1_n1,n2,n3,rnd,holes;
public:
    Supershape(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t n1 = json["n1"].floatValue();
        pm_4 = json["m"].floatValue() / 4.0;
        pneg1_n1 = -1.0 / n1;
        n2 = json["n2"].floatValue();
        n3 = json["n3"].floatValue();
        rnd = json["rnd"].floatValue();
        holes = json["holes"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t theta = pm_4*tx.angle() + M_PI_4;
        num_t st,ct;
        math::sincosg(theta,st,ct);
        num_t t1 = pow(fabs(ct),n2);
        num_t t2 = pow(fabs(st),n3);
        num_t tr = tx.norm2();
        num_t r = (rnd*rng.randNum() + (1.0-rnd)*tr) - holes;
        r *= pow(t1+t2,pneg1_n1) / tr;
        return r * tx;
    }
};

/*
Flower - 2d, from flam3 var51_flower
*/
template <size_t dims>
class Flower: public VariationFrom2D<dims>
{
    num_t petals,holes;
public:
    Flower(const Json& json): VariationFrom2D<dims>(json)
    {
        petals = json["petals"].floatValue();
        holes = json["holes"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t theta = tx.angle();
        num_t r = (rng.randNum() - holes) * cos(petals*theta);
        r /= tx.norm2() + eps_v<num_t>; // flam3 does not add eps
        return r * tx;
    }
};

/*
Conic - 2d, from flam3 var52_conic
*/
template <size_t dims>
class Conic: public VariationFrom2D<dims>
{
    num_t eccen,holes;
public:
    Conic(const Json& json): VariationFrom2D<dims>(json)
    {
        eccen = json["eccen"].floatValue();
        holes = json["holes"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t tr = tx.norm2();
        num_t ct = tx.x() / (tr + eps_v<num_t>); // flam3 does not add eps
        num_t r = (rng.randNum() - holes) * eccen / (tr + tr*eccen*ct);
        return r * tx;
    }
};

/*
Parabola - 2d, from flam3 var53_parabola
*/
template <size_t dims>
class Parabola: public VariationFrom2D<dims>
{
    num_t h,w;
public:
    Parabola(const Json& json): VariationFrom2D<dims>(json)
    {
        h = json["height"].floatValue();
        w = json["width"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t sr,cr;
        math::sincosg(tx.norm2(),sr,cr);
        num_t x = h*sr*sr*rng.randNum();
        num_t y = w*cr*rng.randNum();
        return Point<num_t,2>(x,y);
    }
};

/*
Bipolar - 2d, from flam3 var55_bipolar
*/
template <size_t dims>
class Bipolar: public VariationFrom2D<dims>
{
    num_t shift;
public:
    Bipolar(const Json& json): VariationFrom2D<dims>(json)
    {
        shift = -M_PI_2 * json["shift"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x2y2 = tx.norm2sq();
        num_t t = x2y2 + 1.0;
        num_t x2 = 2.0*tx.x();
        num_t y = 0.5*atan2(2.0*tx.y(),x2y2-1.0) + shift;
        // this should behave the same as flam3
        y -= M_PI * floor(y*M_1_PI + 0.5);
        // not including 0.25 factor of the log that flam3 has
        return Point<num_t,2>(log((t+x2)/(t-x2)),y);
    }
};

/*
Boarders - 2d, from flam3 var56_boarders
Adds a customization to the probability
*/
template <size_t dims>
class Boarders: public VariationFrom2D<dims>
{
    num_t prob;
public:
    Boarders(const Json& json): VariationFrom2D<dims>(json)
    {
        // in flam3 this is 0.75, should this be default?
        prob = json["prob"].floatValue();
        if (prob < 0.0 || prob > 1.0)
            throw std::runtime_error("boarders probability out of range");
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t x,y;
        tx.getXY(x,y);
        num_t rx = rint(x);
        num_t ry = rint(y);
        num_t ox = x - rx;
        num_t oy = y - ry;
        if (rng.randNum() >= prob)
            return Point<num_t,2>(ox*0.5+rx,oy*0.5+ry);
        else
        {
            // in flam3, this is 0.25
            // should it be customizable or use this formula?
            num_t mag = 1.0 - prob;
            if (fabs(ox) >= fabs(oy))
            {
                num_t s = copysign(mag,ox);
                num_t x = ox*0.5 + rx + s;
                num_t y = oy*0.5 + ry + s*oy/ox;
                return Point<num_t,2>(x,y);
            }
            else
            {
                num_t s = copysign(mag,oy);
                num_t x = ox*0.5 + rx + s*ox/oy;
                num_t y = oy*0.5 + ry + s;
                return Point<num_t,2>(x,y);
            }
        }
    }
};

/*
Butterfly - 2d, from flam3 var57_butterfly
*/
template <size_t dims>
class Butterfly: public VariationFrom2D<dims>
{
    // constant from flam3, not used here
    static constexpr num_t flam3_constant = 1.3029400317411197908970256609023;
public:
    Butterfly(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t y2 = 2.0*y;
        num_t r = sqrt(fabs(x*y) / (x*x + y2*y2 + eps_v<num_t>));
        return r * Point<num_t,2>(x,y2);
    }
};

/*
Cell - 2d, from flam3 var58_cell
*/
template <size_t dims>
class Cell: public VariationFrom2D<dims>
{
    num_t size,invsize;
public:
    Cell(const Json& json): VariationFrom2D<dims>(json)
    {
        size = json["size"].floatValue();
        invsize = 1.0/size;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x = floor(tx.x() * invsize);
        num_t y = floor(tx.y() * invsize);
        num_t dx = tx.x() - x*size;
        num_t dy = tx.y() - y*size;
        num_t xs = copysign(2.0,x);
        num_t ys = copysign(2.0,y);
        num_t x2 = x * xs;
        num_t y2 = y * ys;
        // use boolean converted to float
        x2 -= (num_t)(x < 0);
        y2 -= (num_t)(y < 0);
        return Point<num_t,2>(dx+x2*size,-dy-y2*size);
    }
};

/*
CPow - 2d, from flam3 var59_cpow
*/
template <size_t dims>
class CPow: public VariationFrom2D<dims>
{
    num_t va,vc,vd,power;
public:
    CPow(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t r = json["r"].floatValue();
        num_t i = json["i"].floatValue();
        power = json["power"].floatValue();
        va = 2.0*M_PI/power;
        vc = r/power;
        vd = i/power;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        num_t a = tx.angle();
        num_t lnr = 0.5 * log(tx.norm2sq());
        num_t ang = vc*a + vd*lnr + va*floor(power*rng.randNum());
        num_t sa,ca;
        math::sincosg(ang,sa,ca);
        return exp(vc*lnr - vd*a) * Point<num_t,2>(ca,sa);
    }
};

/*
Curve - 2d, from flam3 var60_curve
*/
template <size_t dims>
class Curve: public VariationFrom2D<dims>
{
    num_t invxl,invyl,xamp,yamp;
public:
    Curve(const Json& json): VariationFrom2D<dims>(json)
    {
        xamp = json["xamp"].floatValue();
        yamp = json["yamp"].floatValue();
        num_t xlen = json["xlen"].floatValue();
        num_t ylen = json["ylen"].floatValue();
        xlen *= xlen;
        ylen *= ylen;
        // flam3 uses 1e-20 as the minimum value for the denominator
        invxl = 1.0 / std::max(eps_v<num_t>,xlen);
        invyl = 1.0 / std::max(eps_v<num_t>,ylen);
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t vx = xamp*exp(-y*y*invxl);
        num_t vy = yamp*exp(-x*x*invyl);
        return tx + Point<num_t,2>(vx,vy);
    }
};

/*
E Disc - 2d, from flam3 var61_edisc
*/
template <size_t dims>
class EDisc: public VariationFrom2D<dims>
{
    // constants from flam3, unused
    static constexpr num_t flam3_constant_inv = 0.0864278365005759;
    static constexpr num_t flam3_constant = 11.57034632;
public:
    EDisc(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t tmp = tx.norm2sq() + 1.0;
        num_t tmp2 = 2.0*x;
        num_t xmax = 0.5*(sqrt(tmp+tmp2) + sqrt(tmp-tmp2));
        num_t a1 = log(xmax + sqrt(xmax-1.0));
        num_t a2 = -acos(x/xmax);
        num_t s1,c1;
        math::sincosg(a1,s1,c1);
        num_t s2 = sinh(a2);
        num_t c2 = cosh(a2);
        s1 *= copysign(1.0,-y);
        return Point<num_t,2>(c2*c1,s2*s1);
    }
};

/*
Elliptic - 2d, from flam3 var62_elliptic
*/
template <size_t dims>
struct Elliptic: public VariationFrom2D<dims>
{
    Elliptic(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t tmp = tx.norm2sq() + 1.0;
        num_t x2 = 2.0*x;
        num_t xmax = 0.5*(sqrt(tmp+x2) + sqrt(tmp-x2));
        num_t a = x/xmax;
        num_t b = 1.0 - a*a;
        num_t ssx = xmax - 1.0;
        // is there a way to eliminate these branches?
        b = b < 0.0 ? 0.0 : sqrt(b);
        ssx = ssx < 0.0 ? 0.0 : sqrt(ssx);
        return Point<num_t,2>(atan2(a,b),copysign(1.0,y)*log(xmax+ssx));
    }
};

/*
Escher - 2d, from flam3 var63_escher
*/
template <size_t dims>
class Escher: public VariationFrom2D<dims>
{
    num_t vc,vd;
public:
    Escher(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t beta = json["beta"].floatValue();
        num_t seb,ceb;
        math::sincosg(beta,seb,ceb);
        vc = 0.5*(1.0+ceb);
        vd = 0.5*seb;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t lnr = 0.5*log(tx.norm2sq());
        num_t n = vc*a + vd*lnr;
        num_t sn,cn;
        math::sincosg(n,sn,cn);
        return exp(vc*lnr - vd*a) * Point<num_t,2>(cn,sn);
    }
};

/*
Foci - 2d, from flam3 var64_foci
*/
template <size_t dims>
struct Foci: public VariationFrom2D<dims>
{
    Foci(const Json& json): VariationFrom2D<dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t expx = 0.5*exp(x);
        num_t expnx = 0.25/expx;
        num_t sn,cn;
        math::sincosg(y,sn,cn);
        num_t tmp = 1.0 / (expx + expnx - cn);
        return tmp * Point<num_t,2>(expx-expnx,sn);
    }
};

/*
Lazy Susan - 2d, from flam3 var65_lazysusan
Uses a nonstandard variation weight
*/
template <size_t dims>
class LazySusan: public VariationFrom2D<dims>
{
    num_t px,py,spin,twist,space,lsweight;
public:
    LazySusan(const Json& json): VariationFrom2D<dims>(json)
    {
        px = json["x"].floatValue();
        py = json["y"].floatValue();
        spin = json["spin"].floatValue();
        twist = json["twist"].floatValue();
        space = json["space"].floatValue();
        lsweight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x = tx.x() - px;
        num_t y = tx.y() + py;
        num_t r = hypot(x,y);
        if (r < lsweight)
        {
            num_t a = atan2(y,x) + spin + twist*(lsweight - r);
            num_t sa,ca;
            math::sincosg(a,sa,ca);
            return Point<num_t,2>(r*ca+px,r*sa-py);
        }
        else
        {
            r = 1.0 + space / (r + eps_v<num_t>); // flam3 does not add eps
            return Point<num_t,2>(r*x+px,r*y-py);
        }
    }
};

/*
Loonie - 2d, from flam3 var66_loonie
Uses a nonstandard variation weight
*/
template <size_t dims>
class Loonie: public VariationFrom2D<dims>
{
    num_t loonieweight,loonieweightsq;
public:
    Loonie(const Json& json): VariationFrom2D<dims>(json)
    {
        loonieweight = json["flam3_weight"].floatValue();
        loonieweightsq = loonieweight * loonieweight;
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r2 = tx.norm2sq();
        num_t w2 = loonieweightsq;
        num_t r = loonieweight;
        // flam3 does not add eps
        if (r2 < w2) r *= sqrt(w2/(r2 + eps_v<num_t>) - 1.0);
        return r * tx;
    }
};

/*
O Scope - 2d, from flam3 var69_oscope
*/
template <size_t dims>
class OScope: public VariationFrom2D<dims>
{
    num_t tpf,p_amp,p_damp,sep;
public:
    OScope(const Json& json): VariationFrom2D<dims>(json)
    {
        num_t freq = json["frequency"].floatValue();
        tpf = 2.0*M_PI*freq;
        p_amp = json["amplitude"].floatValue();
        p_damp = json["damping"].floatValue();
        sep = json["separation"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t damp = exp(-fabs(x)*p_damp);
        num_t t = p_amp * damp * cos(tpf*x) + sep;
        num_t yy = copysign(1.0,fabs(y)-t) * y;
        return Point<num_t,2>(x,yy);
    }
};

/*
Popcorn - 2d, from flam3 var17_popcorn and var71_popcorn2
*/
template <size_t dims>
class Popcorn: public VariationFrom2D<dims>
{
    num_t px,py,pc;
public:
    Popcorn(const Json& json): VariationFrom2D<dims>(json)
    {
        px = json["x"].floatValue();
        py = json["y"].floatValue();
        pc = json["c"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t dx = px*sin(tan(tx.y()*pc));
        num_t dy = py*sin(tan(tx.x()*pc));
        return tx + Point<num_t,2>(dx,dy);
    }
};

/*
======= variations by tkoz =======
*/

/*
Spherical P - tkoz, generalized from spherical
*/
template <size_t dims>
class SphericalP: public Variation<dims>
{
    num_t norm;
public:
    SphericalP(const Json& json): Variation<dims>(json)
    {
        norm = json["norm"].floatValue();
        if (norm <= 0.0)
            throw std::runtime_error("norm <= 0");
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.normsum(norm) + eps_v<num_t>);
        return r * tx;
    }
};

/*
Unit Sphere - tkoz
*/
template <size_t dims>
struct UnitSphere: public Variation<dims>
{
    UnitSphere(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2() + eps_v<num_t>);
        return r * tx;
    }
};

/*
Unit Sphere P - tkoz
*/
template <size_t dims>
class UnitSphereP: public Variation<dims>
{
    num_t norm;
public:
    UnitSphereP(const Json& json): Variation<dims>(json)
    {
        norm = json["norm"].floatValue();
        if (norm <= 0.0)
            throw std::runtime_error("norm <= 0");
    }
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm(norm) + eps_v<num_t>);
        return r * tx;
    }
};

/*
Unit Cube - tkoz
*/
template <size_t dims>
struct UnitCube: public Variation<dims>
{
    UnitCube(const Json& json): Variation<dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norminf() + eps_v<num_t>);
        return r * tx;
    }
};

// do not include 2d variations with below 2 dimensions
template <size_t dims>
ENABLE_IF3(dims<2) Variation<dims>::parseVariation2d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put 2d variation factory here
template <size_t dims>
ENABLE_IF3(dims>=2) Variation<dims>::parseVariation2d(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "swirl")
        return new Swirl<dims>(json);
    if (name == "horseshoe")
        return new Horseshoe<dims>(json);
    if (name == "polar")
        return new Polar<dims>(json);
    if (name == "polar2")
        return new Polar2<dims>(json);
    if (name == "handkerchief")
        return new Handkerchief<dims>(json);
    if (name == "heart")
        return new Heart<dims>(json);
    if (name == "disc")
        return new Disc<dims>(json);
    if (name == "disc2")
        return new Disc2<dims>(json);
    if (name == "waves")
        return new Waves<dims>(json);
    if (name == "fan")
        return new Fan<dims>(json);
    if (name == "rings")
        return new Rings<dims>(json);
    if (name == "spiral")
        return new Spiral<dims>(json);
    if (name == "hyperbolic")
        return new Hyperbolic<dims>(json);
    if (name == "diamond")
        return new Diamond<dims>(json);
    if (name == "ex")
        return new Ex<dims>(json);
    if (name == "julia")
        return new Julia<dims>(json);
    if (name == "exponential")
        return new Exponential<dims>(json);
    if (name == "power")
        return new Power<dims>(json);
    if (name == "cosine")
        return new Cosine<dims>(json);
    if (name == "blob")
        return new Blob<dims>(json);
    if (name == "pdj")
        return new PDJ<dims>(json);
    if (name == "cylinder")
        return new Cylinder<dims>(json);
    if (name == "perspective")
        return new Perspective<dims>(json);
    if (name == "julian")
        return new JuliaN<dims>(json);
    if (name == "juliascope")
        return new JuliaScope<dims>(json);
    if (name == "radial_blur")
        return new RadialBlur<dims>(json);
    if (name == "pie")
        return new Pie<dims>(json);
    if (name == "ngon")
        return new NGon<dims>(json);
    if (name == "curl")
        return new Curl<dims>(json);
    if (name == "arch")
        return new Arch<dims>(json);
    if (name == "tangent")
        return new Tangent<dims>(json);
    if (name == "rays")
        return new Rays<dims>(json);
    if (name == "blade")
        return new Blade<dims>(json);
    if (name == "secant")
        return new Secant<dims>(json);
    if (name == "twintrian")
        return new Twintrian<dims>(json);
    if (name == "cross")
        return new Cross<dims>(json);
    if (name == "exp")
        return new Exp<dims>(json);
    if (name == "log")
        return new Log<dims>(json);
    if (name == "sin")
        return new Sin<dims>(json);
    if (name == "cos")
        return new Cos<dims>(json);
    if (name == "tan")
        return new Tan<dims>(json);
    if (name == "sec")
        return new Sec<dims>(json);
    if (name == "csc")
        return new Csc<dims>(json);
    if (name == "cot")
        return new Cot<dims>(json);
    if (name == "sinh")
        return new Sinh<dims>(json);
    if (name == "cosh")
        return new Cosh<dims>(json);
    if (name == "tanh")
        return new Tanh<dims>(json);
    if (name == "sech")
        return new Sech<dims>(json);
    if (name == "csch")
        return new Csch<dims>(json);
    if (name == "coth")
        return new Coth<dims>(json);
    if (name == "auger")
        return new Auger<dims>(json);
    if (name == "flux")
        return new Flux<dims>(json);
    if (name == "mobius")
        return new Mobius<dims>(json);
    if (name == "scry")
        return new Scry<dims>(json);
    if (name == "split")
        return new Split<dims>(json);
    if (name == "stripes")
        return new Stripes<dims>(json);
    if (name == "wedge")
        return new Wedge<dims>(json);
    if (name == "wedge_julia")
        return new WedgeJulia<dims>(json);
    if (name == "wedge_sph")
        return new WedgeSph<dims>(json);
    if (name == "whorl")
        return new Whorl<dims>(json);
    if (name == "supershape")
        return new Supershape<dims>(json);
    if (name == "flower")
        return new Flower<dims>(json);
    if (name == "conic")
        return new Conic<dims>(json);
    if (name == "parabola")
        return new Parabola<dims>(json);
    if (name == "bipolar")
        return new Bipolar<dims>(json);
    if (name == "boarders")
        return new Boarders<dims>(json);
    if (name == "butterfly")
        return new Butterfly<dims>(json);
    if (name == "cell")
        return new Cell<dims>(json);
    if (name == "cpow")
        return new CPow<dims>(json);
    if (name == "curve")
        return new Curve<dims>(json);
    if (name == "edisc")
        return new EDisc<dims>(json);
    if (name == "elliptic")
        return new Elliptic<dims>(json);
    if (name == "escher")
        return new Escher<dims>(json);
    if (name == "foci")
        return new Foci<dims>(json);
    if (name == "lazysusan")
        return new LazySusan<dims>(json);
    if (name == "loonie")
        return new Loonie<dims>(json);
    if (name == "oscope")
        return new OScope<dims>(json);
    if (name == "popcorn")
        return new Popcorn<dims>(json);
    return nullptr;
}

// do not include 3d variations with below 3 dimensions
template <size_t dims>
ENABLE_IF3(dims<3) Variation<dims>::parseVariation3d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put 3d variation factory here
template <size_t dims>
ENABLE_IF3(dims>=3) Variation<dims>::parseVariation3d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put nd variation factory here
template <size_t dims>
Variation<dims> *Variation<dims>::parseVariationNd(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "linear")
        return new Linear<dims>(json);
    if (name == "sinusoidal")
        return new Sinusoidal<dims>(json);
    if (name == "spherical")
        return new Spherical<dims>(json);
    if (name == "bent")
        return new Bent<dims>(json);
    if (name == "rectangles")
        return new Rectangles<dims>(json);
    if (name == "fisheye")
        return new Fisheye<dims>(json);
    if (name == "bubble")
        return new Bubble<dims>(json);
    if (name == "noise")
        return new Noise<dims>(json);
    if (name == "blur")
        return new Blur<dims>(json);
    if (name == "gaussian_blur")
        return new GaussianBlur<dims>(json);
    if (name == "square_noise")
        return new SquareNoise<dims>(json);
    if (name == "separation")
        return new Separation<dims>(json);
    if (name == "splits")
        return new Splits<dims>(json);
    if (name == "pre_blur")
        return new PreBlur<dims>(json);
    if (name == "modulus")
        return new Modulus<dims>(json);
    if (name == "celln")
        return new CellN<dims>(json);
    if (name == "spherical_p")
        return new SphericalP<dims>(json);
    if (name == "unit_sphere")
        return new UnitSphere<dims>(json);
    if (name == "unit_sphere_p")
        return new UnitSphereP<dims>(json);
    if (name == "unit_cube")
        return new UnitCube<dims>(json);
    return nullptr;
}

template <size_t dims>
Variation<dims> *Variation<dims>::parseVariation(const Json& json)
{
    std::string name = json["name"].stringValue();
    Variation<dims> *ret = nullptr;
    ret = parseVariation2d(json);
    if (ret) return ret;
    ret = parseVariation3d(json);
    if (ret) return ret;
    ret = parseVariationNd(json);
    if (ret) return ret;
    throw std::runtime_error("unknown variation: "+name);
}

}

#undef ENABLE_IF
#undef ENABLE_IF2
#undef ENABLE_IF3
