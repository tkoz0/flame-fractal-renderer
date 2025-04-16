/*
Representation of a flame fractal
*/

#pragma once

#include "xform.hpp"

#include "../utils/json.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <vector>

namespace tkoz::flame
{

// flame fractal
template <size_t dims>
class Flame
{
private:

    // number of grid samples in each dimension
    std::array<size_t,dims> size;
    // coordinate bounds in each dimension
    std::array<std::pair<num_t,num_t>,dims> bounds;
    // list of (non final) xforms
    std::vector<XForm<dims>> xforms;
    std::vector<XForm<dims>> xforms_removed;
    // count input xforms because some may be removed in optimization
    size_t xform_ids;
    // final xform
    XForm<dims> final_xform;
    bool has_final_xform;
    // cumulative weights for xform probability selection
    std::vector<num_t> xfcw;
    // color information
    size_t color_dims;
    num_t color_speed;

    // cumulative weights for selecting xform by random number in [0,1)
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

    // make some changes to optimize xform a bit
    void _optimize()
    {
        // remove 0 weight xforms
        auto itr = xforms.begin();
        while (itr != xforms.end())
        {
            if (itr->getWeight() == 0.0)
            {
                xforms_removed.push_back(*itr);
                itr = xforms.erase(itr);
            }
            else
                ++itr;
        }
        if (xforms.empty())
            throw JsonError("Flame(): no xforms remaining after optimization");
        // sort by decreasing weight
        std::sort(xforms.begin(),xforms.end(),
            [](XForm<dims>& a, XForm<dims>& b)
            { return a.getWeight() > b.getWeight(); });
    }

    Flame(){}

public:

    // construct from JSON data
    // throws error if something goes wrong
    Flame(const Json& input)
    {
        try
        {
            JsonArray sizej = input["size"].arrayValue();
            JsonArray boundsj = input["bounds"].arrayValue();
            if (sizej.size() != dims)
                throw JsonError("incorrect size length");
            if (boundsj.size() != dims)
                throw JsonError("incorrect bounds length");
            num_t M = max_rect_v<num_t>;
            for (size_t i = 0; i < dims; ++i)
            {
                try
                {
                    size[i] = sizej[i].floatValue();
                }
                catch (std::exception& e)
                {
                    throw JsonError("error parsing size[" + std::to_string(i)
                        + "]: " + e.what());
                }
                if (size[i] == 0 || size[i] > max_dim)
                    throw JsonError("size[" + std::to_string(i)
                        + "] out of range");
                JsonArray boundj = boundsj[i].arrayValue();
                if (boundj.size() != 2)
                    throw JsonError("bounds[" + std::to_string(i)
                        + "] wrong format");
                num_t lo,hi;
                try
                {
                    lo = boundj[0].floatValue();
                    hi = boundj[1].floatValue();
                }
                catch (std::exception& e)
                {
                    throw JsonError("error parsing bounds[" + std::to_string(i)
                        + "]: " + e.what());
                }
                bounds[i] = std::make_pair(lo,hi);
                if (lo < -M || lo > M || hi < -M || hi > M)
                    throw JsonError("bounds[" + std::to_string(i)
                        + "] out of range");
                if (lo >= hi)
                    throw JsonError("bounds[" + std::to_string(i)
                        + "] low >= high");
            }
        }
        catch (std::exception& e)
        {
            throw JsonError("Flame(): " + std::string(e.what()));
        }
        Json fxf;
        has_final_xform = input.valueAt("final_xform",fxf);
        JsonArray xfs;
        try
        {
            xfs = input["xforms"].arrayValue();
        }
        catch (std::exception& e)
        {
            throw JsonError("Flame(): cannot parse xforms: "
                + std::string(e.what()));
        }
        // color initial
        Json cd,cs;
        try
        {
            if (input.valueAt("color_dimensions",cd))
                color_dims = cd.intValue();
            else
                color_dims = 0;
            if (color_dims > 127)
                throw JsonError("too many color dimensions");
            if (input.valueAt("color_speed",cs))
                color_speed = cs.floatValue();
            else
                color_speed = 0.5;
            if (color_speed < 0.0 || color_speed > 1.0)
                throw std::runtime_error("color speed out of range");
        }
        catch (std::exception& e)
        {
            throw JsonError("Flame(): " + std::string(e.what()));
        }
        // xforms loop
        size_t id = 0;
        for (Json xf : xfs)
        {
            try
            {
                xforms.push_back(XForm<dims>(xf,id,false,color_dims,
                    color_speed));
            }
            catch (std::exception& e)
            {
                throw JsonError("Flame(): error parsing xforms["
                    + std::to_string(id) + "]: " + e.what());
            }
            ++id;
        }
        xform_ids = id;
        if (xforms.empty())
            throw JsonError("Flame(): no xforms");
        if (has_final_xform)
        {
            try
            {
                final_xform = XForm<dims>(fxf,-1,true,color_dims,color_speed);
            }
            catch (std::exception& e)
            {
                throw JsonError("Flame(): error parsing final xform: "
                    + std::string(e.what()));
            }
        }
        _optimize();
        _setupCumulativeWeights();
    }

    inline const XForm<dims>& getRandomXForm(rng_t& rng) const
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

    inline const std::vector<XForm<dims>>& getXForms() const
    {
        return xforms;
    }

    inline const std::vector<XForm<dims>>& getRemovedXForms() const
    {
        return xforms_removed;
    }

    inline size_t getXFormIDCount() const
    {
        return xform_ids;
    }

    inline const XForm<dims>& getFinalXForm() const
    {
        return final_xform;
    }

    inline bool hasFinalXForm() const
    {
        return has_final_xform;
    }

    inline const std::vector<num_t>& getCumulativeWeights() const
    {
        return xfcw;
    }

    inline size_t getColorDims() const
    {
        return color_dims;
    }

    inline num_t getDefaultColorSpeed() const
    {
        return color_speed;
    }
};

} // namespace tkoz::flame
