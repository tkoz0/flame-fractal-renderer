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
public:
    virtual ~Variation(){}
    virtual void calc(IterState<num_t,dims>& state) const = 0;
    static Variation<num_t,dims> *parseVariation(const Json& json);
};

}
