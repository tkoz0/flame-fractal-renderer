/*
Representation of a flame fractal
*/

#pragma once

#include <cstdlib>

#include "xform.hpp"
#include "constants.hpp"

namespace tkoz
{
namespace flame
{

// flame fractal
template <typename num_t, size_t dims>
class Flame
{
private:
    std::array<size_t,dims> size;
    std::array<std::pair<num_t,num_t>,dims> bounds;
    std::vector<XForm<num_t,dims>> xforms;
    XForm<num_t,dims> final_xform;
    bool has_final_xform;
public:
    Flame(){}
    // construct from JSON data
    // throws error if something goes wrong
    Flame(const Json& input)
    {
        if (input["dimensions"].intValue() != dims)
            throw std::runtime_error("flame: dimension mismatch");
        JsonArray sizej = input["size"].arrayValue();
        JsonArray boundsj = input["bounds"].arrayValue();
        if (sizej.size() != dims)
            throw std::runtime_error("flame: incorrect size length");
        if (boundsj.size() != dims)
            throw std::runtime_error("flame: incorrect bounds length");
        num_t M = max_rect<num_t>::value;
        for (size_t i = 0; i < dims; ++i)
        {
            size[i] = sizej[i].floatValue();
            if (size[i] == 0 || size[i] > max_dim)
                throw std::runtime_error("flame: size not in [1,65535]");
            JsonArray boundj = boundsj[i].arrayValue();
            if (boundj.size() != 2)
                throw std::runtime_error("flame: incorrect bound format");
            num_t lo = boundj[0].floatValue();
            num_t hi = boundj[1].floatValue();
            bounds[i] = std::make_pair(lo,hi);
            if (lo < -M || lo > M || hi < -M || hi > M)
                throw std::runtime_error("flame: bound out of range");
            if (lo >= hi)
                throw std::runtime_error("flame: bound low >= bound high");
        }
        // xforms loop
        for (Json xf : input["xforms"].arrayValue())
        {
            if (xf["weight"].floatValue() == 0.0) // ignore 0 weight xforms
                continue;
            xforms.push_back(XForm<num_t,dims>(xf));
        }
        Json xf;
        has_final_xform = input.valueAt("final_xform",xf);
        if (has_final_xform)
            final_xform = XForm<num_t,dims>(input["final_xform"],true);
    }
    void optimize()
    {
        // sort by decreasing weight
        std::sort(xforms.begin(),xforms.end(),
            [](XForm<num_t,dims>& a, XForm<num_t,dims>& b)
            { return a.getWeight() > b.getWeight(); });
        // optimize each xform
        std::for_each(xforms.begin(),xforms.end(),
            [](XForm<num_t,dims>& xf) { xf.optimize(); });
    }
    inline const std::array<size_t,dims>& getSize() const
    {
        return size;
    }
    inline const std::array<std::pair<num_t,num_t>,dims>& getBounds() const
    {
        return bounds;
    }
    inline const std::vector<XForm<num_t,dims>>& getXForms() const
    {
        return xforms;
    }
    inline bool hasFinalXForm() const
    {
        return has_final_xform;
    }
    inline const XForm<num_t,dims>& getFinalXForm() const
    {
        return final_xform;
    }
};

}
}
