#pragma once

#include "variation_base.hpp"

namespace tkoz::flame::vars
{

/*
Linear - generalized from flam3
*/
template <typename num_t, size_t dims>
struct Linear: public Variation<num_t,dims>
{
    Linear(const Json& json): Variation<num_t,dims>(json) {}
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
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
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
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
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norm2sq() + eps<num_t>::value);
        return r * tx;
    }
};

/*
Swirl - 2d, from flam3
*/
template <typename num_t, size_t dims>
struct Swirl: public VariationFrom2D<num_t,dims>
{
    Swirl(const Json& json): VariationFrom2D<num_t,dims>(json) {}
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,2> calc2d(rng_t& rng, const Point<num_t,2>& tx) const
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
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
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
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
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
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
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
    inline Point<num_t,dims> calc(rng_t& rng, const Point<num_t,dims>& tx) const
    {
        (void)rng;
        num_t r = 1.0 / (tx.norminf() + eps<num_t>::value);
        return r * tx;
    }
};

}
