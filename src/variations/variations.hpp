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
        sincosg(r,&sr,&cr);
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
        sincosg(r*a,&sa,&ca);
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
        sincosg(M_PI*r,&sr,&cr);
        return a * Point<num_t,2>(sr,cr);
    }
};

}
