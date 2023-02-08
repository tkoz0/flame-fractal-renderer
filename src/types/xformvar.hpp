/*
Variation information used within an xform
*/

#pragma once

#include <cstdlib>

#include "iter_state.hpp"

namespace tkoz
{
namespace flame
{

template <typename num_t, size_t dims>
struct XFormVar
{
    // function pointer
    std::function<void(IterState<num_t,dims>&,const num_t*)> func;
    size_t index; // index of first variation parameter in varp (XForm class)
    // parameters are taken in order starting from index
    // the varp vector keeps the parameters compactly in memory
};

}
}
