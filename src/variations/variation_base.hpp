#pragma once

#include <cstdlib>
#include <string>

#include "../utils/json_small.hpp"
#include "../utils/math.hpp"
#include "../utils/sfinae.hpp"
#include "../types/point.hpp"
#include "../types/types.hpp"
#include "../types/constants.hpp"

// macros for using SFINAE
#define ENABLE_IF(COND,RET) template <typename RET2 = RET> \
    typename std::enable_if<(COND),RET2>::type
#define ENABLE_IF2(COND) template <typename RET = Variation<num_t,dims>*> \
    typename std::enable_if<(COND),RET>::type

namespace tkoz::flame::vars
{

// variation base class
template <typename num_t, size_t dims>
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
    virtual inline Point<num_t,dims> calc(
        rng_t<num_t>& rng, const Point<num_t,dims>& tx) const = 0;
    inline num_t getWeight() const { return weight; }
    ENABLE_IF2(dims<2) static parseVariation2d(const Json& json);
    ENABLE_IF2(dims>=2) static parseVariation2d(const Json& json);
    ENABLE_IF2(dims<3) static parseVariation3d(const Json& json);
    ENABLE_IF2(dims>=3) static parseVariation3d(const Json& json);
    static Variation<num_t,dims> *parseVariationNd(const Json& json);
    static Variation<num_t,dims> *parseVariation(const Json& json);
};

// class for generalizing 2d variation to higher dimensions
template <typename num_t, size_t dims>
class VariationFrom2D: public Variation<num_t,dims>
{
    static_assert(dims >= 2);
private:
    u32 axis_x,axis_y;
    typedef Point<num_t,dims> _nd;
    typedef Point<num_t,2> _2d;
public:
    VariationFrom2D(const Json& json): Variation<num_t,dims>(json)
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
    inline _nd calc(rng_t<num_t>& rng, const _nd& tx) const
    {
        return calc_h(rng,tx);
    }
    ENABLE_IF(dims==2,_nd) inline calc_h(rng_t<num_t>& rng, const _nd& tx) const
    {
        return calc2d(rng,tx);
    }
    ENABLE_IF(dims>2,_nd) inline calc_h(rng_t<num_t>& rng, const _nd& tx) const
    {
        _2d ret2d = calc2d(rng,_2d(tx[axis_x],tx[axis_y]));
        _nd ret;
        ret[axis_x] = ret2d.x();
        ret[axis_y] = ret2d.y();
        return ret;
    }
    virtual inline _2d calc2d(rng_t<num_t>& rng, const _2d& tx) const = 0;
};

// class for generalizing 3d variation to higher dimensions
template <typename num_t, size_t dims>
class VariationFrom3D: public Variation<num_t,dims>
{
    static_assert(dims >= 3);
private:
    u32 axis_x,axis_y,axis_z;
    typedef Point<num_t,dims> _nd;
    typedef Point<num_t,3> _3d;
public:
    VariationFrom3D(const Json& json): Variation<num_t,dims>(json)
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
    inline _nd calc(rng_t<num_t>& rng, const _nd& tx) const
    {
        return calc_h(rng,tx);
    }
    ENABLE_IF(dims==3,_nd) inline calc_h(rng_t<num_t>& rng, const _nd& tx) const
    {
        return calc2d(rng,tx);
    }
    ENABLE_IF(dims>3,_nd) inline calc_h(rng_t<num_t>& rng, const _nd& tx) const
    {
        _3d ret3d = calc3d(rng,_3d(tx[axis_x],tx[axis_y],tx[axis_z]));
        _nd ret;
        ret[axis_x] = ret3d.x();
        ret[axis_y] = ret3d.y();
        ret[axis_z] = ret3d.z();
        return ret;
    }
    virtual inline _3d calc3d(rng_t<num_t>& rng, const _3d& tx) const = 0;
};

}

#undef ENABLE_IF
#undef ENABLE_IF2
