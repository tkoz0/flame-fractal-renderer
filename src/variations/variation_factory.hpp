#pragma once

#include <cstdlib>

#include "../utils.hpp"
#include "variation_base.hpp"

namespace tkoz::flame::vars
{

template <typename num_t, size_t dims>
Variation<num_t,dims> *Variation<num_t,dims>::parseVariation(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "linear")
        return new Linear<num_t,dims>(json);
    else if (name == "sinusoidal")
        return new Sinusoidal<num_t,dims>(json);
    else if (name == "spherical")
        return new Spherical<num_t,dims>(json);
    else if (name == "swirl")
        return new Swirl<num_t,dims>(json);
    else if (name == "horseshoe")
        return new Horseshoe<num_t,dims>(json);
    else if (name == "polar")
        return new Polar<num_t,dims>(json);
    else if (name == "handkerchief")
        return new Handkerchief<num_t,dims>(json);
    else if (name == "heart")
        return new Heart<num_t,dims>(json);
    else if (name == "disc")
        return new Disc<num_t,dims>(json);
    else if (name == "spiral")
        return new Spiral<num_t,dims>(json);
    else if (name == "hyperbolic")
        return new Hyperbolic<num_t,dims>(json);
    else if (name == "diamond")
        return new Diamond<num_t,dims>(json);
    else if (name == "ex")
        return new Ex<num_t,dims>(json);
    else
        throw std::runtime_error("unknown variation");
}

}
