#pragma once

#include "variation_base.hpp"

namespace tkoz::flame::vars
{

template <typename num_t, size_t dims>
class Spherical: public Variation<num_t,dims>
{
private:
    num_t weight;
public:
    Spherical(const Json& json)
    {
        weight = json["weight"].floatValue();
    }
    void calc(IterState<num_t,dims>& state) const
    {
        num_t r = weight / (state.t.norm2sq() + eps<num_t>::value);
        state.v += r * state.t;
    }
};

}
