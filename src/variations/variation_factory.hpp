#pragma once

#include <cstdlib>

#include "../utils/json_small.hpp"
#include "variation_base.hpp"

// macros for using SFINAE
#define ENABLE_IF2(COND) template <typename RET> \
    typename std::enable_if<(COND),RET>::type

namespace tkoz::flame::vars
{

// do not include 2d variations with below 2 dimensions
template <typename num_t, size_t dims>
ENABLE_IF2(dims<2) Variation<num_t,dims>::parseVariation2d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put 2d variation factory here
template <typename num_t, size_t dims>
ENABLE_IF2(dims>=2) Variation<num_t,dims>::parseVariation2d(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "swirl")
        return new Swirl<num_t,dims>(json);
    if (name == "horseshoe")
        return new Horseshoe<num_t,dims>(json);
    if (name == "polar")
        return new Polar<num_t,dims>(json);
    if (name == "handkerchief")
        return new Handkerchief<num_t,dims>(json);
    if (name == "heart")
        return new Heart<num_t,dims>(json);
    if (name == "disc")
        return new Disc<num_t,dims>(json);
    if (name == "spiral")
        return new Spiral<num_t,dims>(json);
    if (name == "hyperbolic")
        return new Hyperbolic<num_t,dims>(json);
    if (name == "diamond")
        return new Diamond<num_t,dims>(json);
    if (name == "ex")
        return new Ex<num_t,dims>(json);
    if (name == "julia")
        return new Julia<num_t,dims>(json);
    return nullptr;
}

// do not include 3d variations with below 3 dimensions
template <typename num_t, size_t dims>
ENABLE_IF2(dims<3) Variation<num_t,dims>::parseVariation3d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put 3d variation factory here
template <typename num_t, size_t dims>
ENABLE_IF2(dims>=3) Variation<num_t,dims>::parseVariation3d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put nd variation factory here
template <typename num_t, size_t dims>
Variation<num_t,dims> *Variation<num_t,dims>::parseVariationNd(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "linear")
        return new Linear<num_t,dims>(json);
    if (name == "sinusoidal")
        return new Sinusoidal<num_t,dims>(json);
    if (name == "spherical")
        return new Spherical<num_t,dims>(json);
    if (name == "bent")
        return new Bent<num_t,dims>(json);
    if (name == "rectangles")
        return new Rectangles<num_t,dims>(json);
    if (name == "fisheye")
        return new Fisheye<num_t,dims>(json);
    if (name == "bubble")
        return new Bubble<num_t,dims>(json);
    if (name == "noise")
        return new Noise<num_t,dims>(json);
    if (name == "blur")
        return new Blur<num_t,dims>(json);
    if (name == "gaussian_blur")
        return new GaussianBlur<num_t,dims>(json);
    if (name == "square_noise")
        return new SquareNoise<num_t,dims>(json);
    if (name == "spherical_p")
        return new SphericalP<num_t,dims>(json);
    if (name == "unit_sphere")
        return new UnitSphere<num_t,dims>(json);
    if (name == "unit_sphere_p")
        return new UnitSphereP<num_t,dims>(json);
    if (name == "unit_cube")
        return new UnitCube<num_t,dims>(json);
    return nullptr;
}

template <typename num_t, size_t dims>
Variation<num_t,dims> *Variation<num_t,dims>::parseVariation(const Json& json)
{
    std::string name = json["name"].stringValue();
    Variation<num_t,dims> *ret = nullptr;
    ret = parseVariation2d(json);
    if (ret) return ret;
    ret = parseVariation3d(json);
    if (ret) return ret;
    ret = parseVariationNd(json);
    if (ret) return ret;
    throw std::runtime_error("unknown variation: "+name);
}

}

#undef ENABLE_IF2
