#pragma once

#include "variation_base.hpp"

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
    Twintrian(const Json& json): VariationsFrom2D<num_t,dims>(json)
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
    Cross(const Json& json): VariationsFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(
        rng_t<num_t>& rng, const Point<num_t,2>& tx) const
    {
        num_t x,y;
        tx.getXY(x,y);
        num_t s = x*x - y*y;
        num_t r = sqrt(1.0 / (s*s + eps<num_t>::value));
        return r * tx;
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
