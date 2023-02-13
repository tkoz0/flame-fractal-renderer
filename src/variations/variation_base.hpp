#pragma once

#include <cstdlib>
#include <string>

#include "../utils.hpp"
#include "../types.hpp"

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
    virtual Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const = 0;
    inline num_t getWeight() const { return weight; }
    static Variation<num_t,dims> *parseVariation(const Json& json);
};

// class for generalizing 2d variation to higher dimensions
template <typename num_t, size_t dims>
class VariationFrom2D: public Variation<num_t,dims>
{
private:
    u32 axis_x,axis_y;
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
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        if (dims == 2)
            return calc2d(rng,tx);
        else
            return calc2d(rng,Point<num_t,2>(tx[axis_x],tx[axis_y]));
    }
    virtual inline Point<num_t,2> calc2d(
        rng_t& rng, const Point<num_t,2>& tx) const = 0;
};

// class for generalizing 3d variation to higher dimensions
template <typename num_t, size_t dims>
class VariationFrom3D: public Variation<num_t,dims>
{
private:
    u32 axis_x,axis_y,axis_z;
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
    inline Point<num_t,dims> calc(
        rng_t& rng, const Point<num_t,dims>& tx) const
    {
        if (dims == 3)
            return calc3d(rng,tx);
        else
            return calc3d(rng,Point<num_t,3>(tx[axis_x],tx[axis_y],tx[axis_z]));
    }
    virtual inline Point<num_t,3> calc3d(
        rng_t& rng, const Point<num_t,3>& tx) const = 0;
};

}
