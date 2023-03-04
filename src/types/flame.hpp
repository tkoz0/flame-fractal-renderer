/*
Representation of a flame fractal
*/

#pragma once

#include <algorithm>
#include <array>
#include <cstdlib>
#include <vector>

// forward declaration
namespace tkoz::flame
{
template <typename num_t, size_t dims> class Flame;
}

#include "../utils/json_small.hpp"
#include "xform.hpp"

namespace tkoz::flame
{

// flame fractal
template <typename num_t, size_t dims>
class Flame
{
private:
    // number of grid samples in each dimension
    std::array<size_t,dims> size;
    // coordinate bounds in each dimension
    std::array<std::pair<num_t,num_t>,dims> bounds;
    // list of (non final) xforms
    std::vector<XForm<num_t,dims>> xforms;
    size_t xform_ids;
    // final xform
    XForm<num_t,dims> final_xform;
    bool has_final_xform;
    // cumulative weights for xform probability selection
    std::vector<num_t> xfcw;
    void _setupCumulativeWeights()
    {
        xfcw = std::vector<num_t>(xforms.size());
        // compute weight sum for normalizing
        num_t normdiv = 0.0;
        for (size_t i = 0; i < xforms.size(); ++i)
            normdiv += xforms[i].getWeight();
        // store normalized cumulative weights
        num_t wsum = 0.0;
        for (size_t i = 0; i < xforms.size(); ++i)
        {
            wsum += xforms[i].getWeight() / normdiv;
            xfcw[i] = wsum;
        }
        // correct rounding error
        xfcw.back() = 1.0;
    }
    void _optimize()
    {
        // remove 0 weight xforms
        auto itr = xforms.begin();
        while (itr != xforms.end())
        {
            if (itr->getWeight() == 0.0)
                itr = xforms.erase(itr);
            else
                ++itr;
        }
        // sort by decreasing weight
        std::sort(xforms.begin(),xforms.end(),
            [](XForm<num_t,dims>& a, XForm<num_t,dims>& b)
            { return a.getWeight() > b.getWeight(); });
    }
    Flame(){}
public:
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
                throw std::runtime_error("flame: size out of range");
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
        size_t id = 0;
        for (Json xf : input["xforms"].arrayValue())
        {
            xforms.push_back(XForm<num_t,dims>(xf,id,false));
            ++id;
        }
        xform_ids = id;
        if (xforms.empty())
            throw std::runtime_error("flame: must have >= 1 xform");
        Json xf;
        has_final_xform = input.valueAt("final_xform",xf);
        if (has_final_xform)
            final_xform = XForm<num_t,dims>(input["final_xform"],-1,true);
        _optimize();
        _setupCumulativeWeights();
    }
    inline const XForm<num_t,dims>& getRandomXForm(rng_t<num_t>& rng) const
    {
        size_t i = 0;
        num_t r = rng.randNum();
        while (xfcw[i] < r)
            ++i;
        return xforms[i];
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
    inline size_t getXFormIDCount() const
    {
        return xform_ids;
    }
};

}
