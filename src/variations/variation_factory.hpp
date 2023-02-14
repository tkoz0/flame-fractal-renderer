#pragma once

#include <cstdlib>

#include "../utils/json_small.hpp"
#include "variation_base.hpp"

namespace tkoz::flame::vars
{

template <typename num_t, size_t dims>
Variation<num_t,dims> *Variation<num_t,dims>::parseVariation(const Json& json)
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
    if (name == "spherical_p")
        return new SphericalP<num_t,dims>(json);
    if (name == "unit_sphere")
        return new UnitSphere<num_t,dims>(json);
    if (name == "unit_sphere_p")
        return new UnitSphereP<num_t,dims>(json);
    if (name == "unit_cube")
        return new UnitCube<num_t,dims>(json);
    throw std::runtime_error("unknown variation: "+name);
}

}
