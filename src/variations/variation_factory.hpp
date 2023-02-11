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
    else
        throw std::runtime_error("unknown variation");
}

}
