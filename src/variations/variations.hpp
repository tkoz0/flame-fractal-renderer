#pragma once

#include "variation_base.hpp"
#include "../utils/flame.hpp"

namespace tkoz::flame::vars
{

/*
======= N dimensional variations from flam3 =======
*/

/*
Linear - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Linear: public Variation<num_t,dims>
{
    Linear(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        return tx;
    }
};

/*
Sinusoidal - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Sinusoidal: public Variation<num_t,dims>
{
    Sinusoidal(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        return tx.map([](num_t x){ return sin(x); });
    }
};

/*
Spherical - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Spherical: public Variation<num_t,dims>
{
    Spherical(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2sq() + eps<num_t>::value);
        return r * tx;
    }
};

/*
Bent - generalized from bent2 in flam3
*/
template <typename num_t, size_t dims>
class Bent: public Variation<num_t,dims>
{
    Point<num_t,dims> scales_neg;
    Point<num_t,dims> scales_pos;
public:
    Bent(const Json& json): Variation<num_t,dims>(json)
    {
        scales_neg = Point<num_t,dims>(json["scales_neg"]);
        scales_pos = Point<num_t,dims>(json["scales_pos"]);
    }
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        Point<num_t,dims> ret = tx;
        for (size_t i = 0; i < dims; ++i)
        {
            if (ret[i] < 0.0)
                ret[i] *= scales_neg[i];
            else
                ret[i] *= scales_pos[i];
        }
        return ret;
    }
};

/*
Rectangles - generalized from flam3
*/
template <typename num_t, size_t dims>
class Rectangles: public Variation<num_t,dims>
{
    Point<num_t,dims> params;
public:
    Rectangles(const Json& json): Variation<num_t,dims>(json)
    {
        params = Point<num_t,dims>(json["params"]);
    }
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
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
Fisheye - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Fisheye: public Variation<num_t,dims>
{
    Fisheye(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        // TODO support changing the +1
        num_t r = 1.0 / (tx.norm2() + 1.0);
        return r * tx;
    }
};

/*
Bubble - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Bubble: public Variation<num_t,dims>
{
    Bubble(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        // TODO support changing the +4
        num_t r = 1.0 / (tx.norm2sq() + 4.0);
        return r * tx;
    }
};

/*
Noise - generalized from flamm3
*/
template <typename num_t, size_t dims>
struct Noise: public Variation<num_t,dims>
{
    Noise(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        num_t r = rng.randNum();
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return r * multComponents(tx,dir);
    }
};

/*
Blur - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Blur: public Variation<num_t,dims>
{
    Blur(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        num_t r = rng.randNum();
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return r * dir;
    }
};

/*
Gaussian Blur - generalized from flam3
*/
template <typename num_t, size_t dims>
struct GaussianBlur: public Variation<num_t,dims>
{
    GaussianBlur(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        num_t r = rng.randGaussian();
        Point<num_t,dims> dir = rng.template randDirection<dims>();
        return r * dir;
    }
};

/*
Square Noise - generalized from square in flam3
*/
template <typename num_t, size_t dims>
struct SquareNoise: public Variation<num_t,dims>
{
    SquareNoise(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)tx;
        return rng.template randPoint2<dims>();
    }
};

/*
Separation - generalized from flam3
*/
template <typename num_t, size_t dims>
class Separation: public Variation<num_t,dims>
{
    Point<num_t,dims> sep2,inside;
public:
    Separation(const Json& json): Variation<num_t,dims>(json)
    {
        sep2 = Point<num_t,dims>(json["params"]);
        for (size_t i = 0; i < dims; ++i)
            sep2[i] *= sep2[i];
        inside = Point<num_t,dims>(json["inside"]);
    }
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
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
Splits - generalized from flam3
*/
template <typename num_t, size_t dims>
class Splits: public Variation<num_t,dims>
{
    Point<num_t,dims> params;
public:
    Splits(const Json& json): Variation<num_t,dims>(json)
    {
        params = Point<num_t,dims>(json["params"]);
    }
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
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
======= 2 dimensional variations from flam3 =======
*/

/*
Swirl - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Swirl: public VariationFrom2D<num_t,dims>
{
    Swirl(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = tx.norm2sq();
        num_t x,y;
        tx.getXY(x,y);
        num_t sr,cr;
        sincosg(r,sr,cr);
        return Point<num_t,2>(x*sr-y*cr,x*cr+y*sr);
    }
};

/*
Horseshoe - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Horseshoe: public VariationFrom2D<num_t,dims>
{
    Horseshoe(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2() + eps<num_t>::value);
        num_t x,y;
        tx.getXY(x,y);
        return r * Point<num_t,2>((x-y)*(x+y),2.0*x*y);
    }
};

/*
Polar - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Polar: public VariationFrom2D<num_t,dims>
{
    Polar(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        return Point<num_t,2>(a*M_1_PI,r-1.0);
    }
};

/*
Polar2 - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Polar2: public VariationFrom2D<num_t,dims>
{
    Polar2(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        return Point<num_t,2>(tx.angle(),log(tx.norm2sq()));
    }
};

/*
Handkerchief - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Handkerchief: public VariationFrom2D<num_t,dims>
{
    Handkerchief(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        return r * Point<num_t,2>(sin(a+r),cos(a-r));
    }
};

/*
Heart - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Heart: public VariationFrom2D<num_t,dims>
{
    Heart(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        num_t sa,ca;
        sincosg(r*a,sa,ca);
        return r * Point<num_t,2>(sa,-ca);
    }
};

/*
Disc - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Disc: public VariationFrom2D<num_t,dims>
{
    Disc(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = tx.angle();
        num_t r = tx.norm2();
        num_t sr,cr;
        sincosg(M_PI*r,sr,cr);
        return a * Point<num_t,2>(sr,cr);
    }
};

/*
Spiral - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Spiral: public VariationFrom2D<num_t,dims>
{
    Spiral(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        num_t sr,cr;
        sincosg(r,sr,cr);
        num_t r1 = 1.0 / (r + eps<num_t>::value);
        return r1 * Point<num_t,2>(ca+sr,sa-cr);
    }
};

/*
Hyperbolic - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Hyperbolic: public VariationFrom2D<num_t,dims>
{
    Hyperbolic(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        return Point<num_t,2>(sa/r,ca*(r+eps<num_t>::value));
    }
};

/*
Diamond - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Diamond: public VariationFrom2D<num_t,dims>
{
    Diamond(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        num_t sr,cr;
        sincosg(r,sr,cr);
        return Point<num_t,2>(sa*cr,ca*sr);
    }
};

/*
Ex - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Ex: public VariationFrom2D<num_t,dims>
{
    Ex(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
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
Julia - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Julia: public VariationFrom2D<num_t,dims>
{
    Julia(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        num_t a = 0.5*tx.angle() + rng.randBool()*M_PI; // + 0 or PI
        num_t sa,ca;
        sincosg(a,sa,ca);
        return tx.norm2() * Point<num_t,2>(ca,sa);
    }
};

/*
Exponential - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Exponential: public VariationFrom2D<num_t,dims>
{
    Exponential(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t dx = exp(tx.x()-1.0);
        num_t sdy,cdy;
        sincosg(M_PI*tx.y(),sdy,cdy);
        return dx * Point<num_t,2>(cdy,sdy);
    }
};

/*
Power - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Power: public VariationFrom2D<num_t,dims>
{
    Power(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        return pow(r,sa) * Point<num_t,2>(ca,sa);
    }
};

/*
Cosine - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Cosine: public VariationFrom2D<num_t,dims>
{
    Cosine(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t sa,ca;
        sincosg(x*M_PI,sa,ca);
        return Point<num_t,2>(ca*cosh(y),-sa*sinh(y));
    }
};

/*
Blob - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Blob: public VariationFrom2D<num_t,dims>
{
    num_t mid,amp,waves;
public:
    Blob(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        num_t low = json["low"].floatValue();
        num_t high = json["high"].floatValue();
        mid = (high+low)/2.0;
        amp = (high-low)/2.0;
        waves = json["waves"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r,sa,ca;
        tx.getRadiusSinCos(r,sa,ca);
        num_t a = tx.angle();
        r *= mid + amp*sin(waves*a);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
PDJ - 2d, from flam3
*/
template <typename num_t, size_t dims>
class PDJ: public VariationFrom2D<num_t,dims>
{
    num_t a,b,c,d;
public:
    PDJ(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        a = json["a"].floatValue();
        b = json["b"].floatValue();
        c = json["c"].floatValue();
        d = json["d"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
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
Cylinder - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Cylinder: public VariationFrom2D<num_t,dims>
{
    Cylinder(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        return Point<num_t,2>(sin(tx.x()),tx.y());
    }
};

/*
Perspective - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Perspective: public VariationFrom2D<num_t,dims>
{
    num_t dist,vsin,vfcos;
public:
    Perspective(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        dist = json["distance"].floatValue();
        num_t angle = json["angle"].floatValue();
        vsin = sin(angle);
        vfcos = dist*cos(angle);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t t = 1.0 / (dist - tx.y()*vsin);
        return t * Point<num_t,2>(dist*tx.x(),vfcos*tx.y());
    }
};

/*
Julia N - 2d, from flam3
*/
template <typename num_t, size_t dims>
class JuliaN: public VariationFrom2D<num_t,dims>
{
    num_t abspower,invpower,cn;
public:
    JuliaN(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        num_t power = json["power"].floatValue();
        num_t dist = json["dist"].floatValue();
        abspower = fabs(power);
        invpower = 1.0/power;
        cn = dist/(2.0*power);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        i32 t = trunc(abspower*rng.randNum());
        num_t a = (tx.angle() + (2.0*M_PI)*t) * invpower;
        num_t r = pow(tx.norm2sq(),cn);
        num_t sa,ca;
        sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Julia Scope - 2d, from flam3
*/
template <typename num_t, size_t dims>
class JuliaScope: public VariationFrom2D<num_t,dims>
{
    num_t abspower,invpower,cn;
public:
    JuliaScope(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        num_t power = json["power"].floatValue();
        num_t dist = json["dist"].floatValue();
        abspower = fabs(power);
        invpower = 1.0/power;
        cn = dist/(2.0*power);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        i32 t = trunc(abspower*rng.randNum());
        num_t dir = rng.template randDirection<1>().x();
        num_t a = ((2.0*M_PI)*t + dir*tx.angle()) * invpower;
        num_t r = pow(tx.norm2sq(),cn);
        num_t sa,ca;
        sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Radial Blur - 2d, from flam3
*/
template <typename num_t, size_t dims>
class RadialBlur: public VariationFrom2D<num_t,dims>
{
    num_t spin,zoom,flam3weight;
public:
    RadialBlur(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        num_t angle = json["angle"].floatValue();
        sincosg(angle*M_PI_2,spin,zoom);
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        num_t g = flam3weight * rng.randGaussian();
        num_t ra = tx.norm2();
        num_t a = tx.angle() + spin*g;
        num_t sa,ca;
        sincosg(a,sa,ca);
        num_t rz = zoom*g - 1.0;
        return ra*Point<num_t,2>(ca,sa) + rz*tx;
    }
};

/*
Pie - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Pie: public VariationFrom2D<num_t,dims>
{
    num_t slices,rotation,thickness,invslices2pi;
public:
    Pie(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        slices = json["slices"].floatValue();
        rotation = json["rotation"].floatValue();
        thickness = json["thickness"].floatValue();
        invslices2pi = (2.0*M_PI)/slices;
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)tx;
        i32 sl = (i32)(rng.randNum()*slices + 0.5);
        num_t a = rotation + (sl + rng.randNum()*thickness)*invslices2pi;
        num_t r = rng.randNum();
        num_t sa,ca;
        sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
NGon - 2d, from flam3
*/
template <typename num_t, size_t dims>
class NGon: public VariationFrom2D<num_t,dims>
{
    num_t powerval,angle,corners,circle,invangle;
public:
    NGon(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        num_t sides = json["sides"].floatValue();
        powerval = json["power"].floatValue()/2.0;
        angle = (2.0*M_PI)/sides;
        corners = json["corners"].floatValue();
        circle = json["circle"].floatValue();
        invangle = sides/(2.0*M_PI);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = pow(tx.norm2sq(),powerval);
        num_t theta = tx.angle();
        num_t phi = theta - angle*floor(theta*invangle);
        static const num_t mult[2] = {0.0,1.0};
        phi -= mult[phi > angle*0.5]*angle;
        num_t amp = corners*(1.0/(cos(phi)+eps<num_t>::value)-1.0) + circle;
        amp /= r + eps<num_t>::value;
        return amp * tx;
    }
};

/*
Curl - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Curl: public VariationFrom2D<num_t,dims>
{
    num_t c1,c2;
public:
    Curl(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        c1 = json["c1"].floatValue();
        c2 = json["c2"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t re = 1.0 + c1*x + c2*(x*x - y*y);
        num_t im = c1*y + 2.0*c2*x*y;
        num_t r = 1.0 / (re*re + im*im + eps<num_t>::value);
        return r * Point<num_t,2>(x*re+y*im,y*re-x*im);
    }
};

/*
Arch - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Arch: public VariationFrom2D<num_t,dims>
{
    num_t flam3weight;
public:
    Arch(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)tx;
        num_t a = flam3weight * rng.randNum() * M_PI;
        num_t sa,ca;
        sincosg(a,sa,ca);
        return flam3weight * Point<num_t,2>(sa,sa*sa/ca);
    }
};

/*
Tangent - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Tangent: public VariationFrom2D<num_t,dims>
{
    Tangent(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        return Point<num_t,2>(sin(x)/cos(y),tan(y));
    }
};

/*
Rays - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Rays: public VariationFrom2D<num_t,dims>
{
    num_t flam3weight;
public:
    Rays(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t a = flam3weight * rng.randNum() * M_PI;
        num_t r = flam3weight / (tx.norm2sq() + eps<num_t>::value);
        num_t tr = tan(a) * r;
        return tr * Point<num_t,2>(cos(tx.x()),sin(tx.y()));
    }
};

/*
Blade - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Blade: public VariationFrom2D<num_t,dims>
{
    num_t flam3weight;
public:
    Blade(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        num_t r = rng.randNum() * flam3weight * tx.norm2();
        num_t sr,cr;
        sincosg(r,sr,cr);
        return tx.x() * Point<num_t,2>(cr+sr,cr-sr);
    }
};

/*
Secant - 2d, based on secant2 from flam3
*/
template <typename num_t, size_t dims>
class Secant: public VariationFrom2D<num_t,dims>
{
    num_t flam3weight;
public:
    Secant(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t cr = cos(flam3weight*tx.norm2());
        num_t icr = 1.0/cr;
        static const num_t sign[2] = {-1.0,1.0};
        return Point<num_t,2>(tx.x(),icr+sign[cr<0]);
    }
};

/*
Twintrian - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Twintrian: public VariationFrom2D<num_t,dims>
{
    num_t flam3weight;
public:
    Twintrian(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        flam3weight = json["flam3_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        num_t r = rng.randNum() * flam3weight * tx.norm2();
        num_t sr,cr;
        sincosg(r,sr,cr);
        num_t diff = log10(sr*sr) + cr;
        if (unlikely(bad_value(diff)))
            diff = -30.0;
        return tx.x() * Point<num_t,2>(diff,diff-sr*M_PI);
    }
};

/*
Cross - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Cross: public VariationFrom2D<num_t,dims>
{
    Cross(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s = x*x - y*y;
        num_t r = sqrt(1.0 / (s*s + eps<num_t>::value));
        return r * tx;
    }
};

/*
Exp - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Exp: public VariationFrom2D<num_t,dims>
{
    Exp(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t e = exp(tx.x());
        num_t es,ec;
        sincosg(tx.y(),es,ec);
        return e * Point<num_t,2>(ec,es);
    }
};

/*
Log - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Log: public VariationFrom2D<num_t,dims>
{
    Log(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        return Point<num_t,2>(log(tx.norm2sq()),tx.angle());
    }
};

/*
Sin - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Sin: public VariationFrom2D<num_t,dims>
{
    Sin(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(x,s,c);
        num_t sh = sinh(y);
        num_t ch = cosh(y);
        return Point<num_t,2>(s*ch,c*sh);
    }
};

/*
Cos - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Cos: public VariationFrom2D<num_t,dims>
{
    Cos(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(x,s,c);
        num_t ch = cosh(y);
        num_t sh = sinh(y);
        return Point<num_t,2>(c*ch,-s*sh);
    }
};

/*
Tan - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Tan: public VariationFrom2D<num_t,dims>
{
    Tan(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(2.0*x,s,c);
        num_t sh = sinh(2.0*y);
        num_t ch = cosh(2.0*y);
        return Point<num_t,2>(s,sh) / (c+ch);
    }
};

/*
Sec - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Sec: public VariationFrom2D<num_t,dims>
{
    Sec(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(x,s,c);
        num_t sh = sinh(y);
        num_t ch = cosh(y);
        return Point<num_t,2>(c*ch,s*sh) / (cos(2.0*x)+cosh(2.0*y));
    }
};

/*
Csc - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Csc: public VariationFrom2D<num_t,dims>
{
    Csc(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(x,s,c);
        num_t sh = sinh(y);
        num_t ch = cosh(y);
        return Point<num_t,2>(s*ch,-c*sh) / (cosh(2.0*y)-cos(2.0*x));
    }
};

/*
Cot - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Cot: public VariationFrom2D<num_t,dims>
{
    Cot(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(2.0*x,s,c);
        num_t sh = sinh(2.0*y);
        num_t ch = cosh(2.0*y);
        return Point<num_t,2>(s,-sh) / (ch-c);
    }
};

/*
Sinh - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Sinh: public VariationFrom2D<num_t,dims>
{
    Sinh(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        return Point<num_t,2>(sh*c,ch*s);
    }
};

/*
Cosh - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Cosh: public VariationFrom2D<num_t,dims>
{
    Cosh(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        return Point<num_t,2>(ch*c,sh*s);
    }
};

/*
Tanh - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Tanh: public VariationFrom2D<num_t,dims>
{
    Tanh(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(2.0*y,s,c);
        num_t sh = sinh(2.0*x);
        num_t ch = cosh(2.0*x);
        return Point<num_t,2>(sh,s) / (c+ch);
    }
};

/*
Sech - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Sech: public VariationFrom2D<num_t,dims>
{
    Sech(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        return Point<num_t,2>(c*ch,-s*sh) / (cos(2.0*y)+cosh(2.0*x));
    }
};

/*
Csch - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Csch: public VariationFrom2D<num_t,dims>
{
    Csch(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(y,s,c);
        num_t sh = sinh(x);
        num_t ch = cosh(x);
        return Point<num_t,2>(sh*c,-ch*s) / (cosh(2.0*x)-cos(2.0*y));
    }
};

/*
Coth - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Coth: public VariationFrom2D<num_t,dims>
{
    Coth(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t s,c;
        sincosg(2.0*y,s,c);
        num_t sh = sinh(2.0*x);
        num_t ch = cosh(2.0*x);
        return Point<num_t,2>(sh,s) / (ch-c);
    }
};

/*
Auger - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Auger: public VariationFrom2D<num_t,dims>
{
    num_t freq,augerweight,scale,sym;
public:
    Auger(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        freq = json["freq"].floatValue();
        augerweight = json["auger_weight"].floatValue();
        scale = json["scale"].floatValue() / 2.0;
        sym = json["sym"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
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
Flux - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Flux: public VariationFrom2D<num_t,dims>
{
    num_t spread,fluxweight;
public:
    Flux(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        spread = 2.0 + json["spread"].floatValue();
        fluxweight = json["flux_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
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
        sincosg(avga,c,s);
        return avgr * Point<num_t,2>(c,s);
    }
};

/*
Mobius - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Mobius: public VariationFrom2D<num_t,dims>
{
    Point<num_t,2> a,b,c,d;
public:
    Mobius(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        a = Point<num_t,2>(json["a"]);
        b = Point<num_t,2>(json["b"]);
        c = Point<num_t,2>(json["c"]);
        d = Point<num_t,2>(json["d"]);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t x,y;
        tx.getXY(x,y);
        num_t re_u = a.x()*x - a.y()*y + b.x();
        num_t im_u = a.x()*y + a.y()*x + b.y();
        num_t re_v = c.x()*x - c.y()*y + d.x();
        num_t im_v = c.x()*y + c.y()*x + d.y();
        num_t rad = 1.0 / (re_v*re_v + im_v*im_v + eps<num_t>::value);
        return rad * Point<num_t,2>(re_u*re_v+im_u*im_v,im_u*re_v-re_u*im_v);
    }
};

/*
Scry - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Scry: public VariationFrom2D<num_t,dims>
{
    num_t scryweight;
public:
    Scry(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        scryweight = json["scry_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t t = tx.norm2sq();
        num_t r = 1.0 / (sqrt(t) * (t + 1.0/(scryweight + eps<num_t>::value)));
        return r * tx;
    }
};

/*
Split - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Split: public VariationFrom2D<num_t,dims>
{
    num_t sizex,sizey;
public:
    Split(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        sizex = json["xsize"].floatValue() * M_PI;
        sizey = json["ysize"].floatValue() * M_PI;
    }
    inline Point<num_t,dims> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t xs = copysign(1.0,cos(tx.x()*sizex));
        num_t ys = copysign(1.0,cos(tx.y()*sizey));
        return Point<num_t,2>(xs*tx.y(),ys*tx.x());
    }
};

/*
Stripes - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Stripes: public VariationFrom2D<num_t,dims>
{
    num_t space,warp;
public:
    Stripes(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        space = 1.0 - json["space"].floatValue();
        warp = json["warp"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t rx = floor(tx.x() + 0.5);
        num_t ox = tx.x() - rx;
        return Point<num_t,2>(ox*space+rx,tx.y()+ox*ox*warp);
    }
};

/*
Wedge - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Wedge: public VariationFrom2D<num_t,dims>
{
    num_t swirl,count,angle,hole;
public:
    Wedge(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        swirl = json["swirl"].floatValue();
        count = json["count"].floatValue();
        angle = json["angle"].floatValue();
        hole = json["hole"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = tx.norm2();
        num_t a = tx.angle() + swirl*r;
        num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
        num_t cf = 1.0 - angle*swirl*(M_1_PI*0.5);
        a = a*cf + c*angle;
        num_t sa,ca;
        sincosg(a,sa,ca);
        return (r+hole) * Point<num_t,2>(ca,sa);
    }
};

/*
Wedge Julia - 2d, from flam3
*/
template <typename num_t, size_t dims>
class WedgeJulia: public VariationFrom2D<num_t,dims>
{
    num_t cn,abspower,invpower,count,angle,cf;
public:
    WedgeJulia(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        angle = json["angle"].floatValue();
        count = json["count"].floatValue();
        num_t power = json["power"].floatValue();
        invpower = 1.0/power;
        num_t dist = json["dist"].floatValue();
        cn = dist/(2.0*power);
        abspower = fabs(power);
        cf = 1.0 - angle*count*(M_1_PI*0.5);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        num_t r = pow(tx.norm2sq(),cn);
        i32 tr = (i32)(abspower * rng.randNum());
        num_t a = (tx.angle() + (2.0*M_PI)*tr) * invpower;
        num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
        num_t sa,ca;
        a = a*cf + c*angle;
        sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
Wedge Sph - 2d, from flam3
*/
template <typename num_t, size_t dims>
class WedgeSph: public VariationFrom2D<num_t,dims>
{
    num_t swirl,count,cf,angle,hole;
public:
    WedgeSph(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        angle = json["angle"].floatValue();
        count = json["count"].floatValue();
        swirl = json["swirl"].floatValue();
        hole = json["hole"].floatValue();
        cf = 1.0 - angle*count*(M_1_PI*0.5);
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2() + eps<num_t>::value);
        num_t a = tx.angle() + swirl*r;
        num_t c = floor((count*a + M_PI) * (M_1_PI*0.5));
        num_t sa,ca;
        a = a*cf + c*angle;
        sincosg(a,sa,ca);
        return (r+hole) * Point<num_t,2>(ca,sa);
    }
};

/*
Whorl - 2d, from flam3
*/
template <typename num_t, size_t dims>
class Whorl: public VariationFrom2D<num_t,dims>
{
    num_t choice[2],whorlweight;
public:
    Whorl(const Json& json): VariationFrom2D<num_t,dims>(json)
    {
        choice[0] = json["inside"].floatValue();
        choice[1] = json["outside"].floatValue();
        whorlweight = json["whorl_weight"].floatValue();
    }
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        (void)rng;
        num_t r = tx.norm2();
        num_t a = tx.angle();
        a += choice[r >= whorlweight] / (whorlweight - r);
        num_t sa,ca;
        sincosg(a,sa,ca);
        return r * Point<num_t,2>(ca,sa);
    }
};

/*
======= variations by tkoz =======
*/

/*
Spherical P - tkoz, generalized from spherical
*/
template <typename num_t, size_t dims>
class SphericalP: public Variation<num_t,dims>
{
    num_t norm;
public:
    SphericalP(const Json& json): Variation<num_t,dims>(json)
    {
        norm = json["norm"].floatValue();
        if (norm <= 0.0)
            throw std::runtime_error("norm <= 0");
    }
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.normsum(norm) + eps<num_t>::value);
        return r * tx;
    }
};

/*
Unit Sphere - tkoz
*/
template <typename num_t, size_t dims>
struct UnitSphere: public Variation<num_t,dims>
{
    UnitSphere(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2() + eps<num_t>::value);
        return r * tx;
    }
};

/*
Unit Sphere P - tkoz
*/
template <typename num_t, size_t dims>
class UnitSphereP: public Variation<num_t,dims>
{
    num_t norm;
public:
    UnitSphereP(const Json& json): Variation<num_t,dims>(json)
    {
        norm = json["norm"].floatValue();
        if (norm <= 0.0)
            throw std::runtime_error("norm <= 0");
    }
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm(norm) + eps<num_t>::value);
        return r * tx;
    }
};

/*
Unit Cube - tkoz
*/
template <typename num_t, size_t dims>
struct UnitCube: public Variation<num_t,dims>
{
    UnitCube(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norminf() + eps<num_t>::value);
        return r * tx;
    }
};

}
