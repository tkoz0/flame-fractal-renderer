#pragma once

#include "variation_base.hpp"

namespace tkoz::flame::vars
{

template <typename num_t, size_t dims>
class Linear: public Variation<num_t,dims>
{
private:
    num_t weight;
public:
    Linear(const Json& json)
    {
        weight = json["weight"].floatValue();
    }
    void calc(IterState<num_t,dims>& state) const
    {
        state.v += weight * state.t;
    }
};

}
