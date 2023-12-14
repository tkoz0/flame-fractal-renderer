/*
Representation of an xform
This class is coupled with Flame and not designed to be used independently
*/

#pragma once

#include <cstdlib>
#include <vector>

// forward declaration
namespace tkoz::flame
{
template <typename num_t, size_t dims> class XForm;
}

#include "point.hpp"
#include "affine.hpp"
#include "../variations.hpp"

namespace tkoz::flame
{

// xform (including final xform)
template <typename num_t, size_t dims>
class XForm
{
private:
    typedef Point<num_t,dims> point_t;
    size_t id;
    num_t weight; // xform probability weight, not applicable for final xform
    typedef vars::Variation<num_t,dims> var_t;
    std::vector<std::shared_ptr<var_t>> vars;
    Affine<num_t,dims> pre; // pre affine transformation
    Affine<num_t,dims> post; // post affine transformation
    bool has_pre,has_post;
    // color vector, length info from flame object
    std::vector<num_t> color;
    // optimize xform
    void _optimize()
    {
        // remove 0 weight variations
        auto itr = vars.begin();
        while (itr != vars.end())
        {
            if (fabs((*itr)->getWeight()) == 0.0)
                itr = vars.erase(itr);
            else
                ++itr;
        }
    }
public:
    // need default constructor since final xform may not be used
    XForm(){}
    // construct from JSON data
    // throws error if something goes wrong
    XForm(const Json& input, size_t id, bool is_final, size_t color_dims)
    {
        this->id = id;
        if (!is_final)
            weight = input["weight"].floatValue();
        else
            weight = 1.0; // unused
        if (weight < 0.0)
            throw std::runtime_error("weights must be nonnegative");
        Json affine;
        has_pre = input.valueAt("pre_affine",affine);
        if (has_pre)
            pre = Affine<num_t,dims>(affine);
        else
            pre = Affine<num_t,dims>();
        has_post = input.valueAt("post_affine",affine);
        if (has_post)
            post = Affine<num_t,dims>(affine);
        else
            post = Affine<num_t,dims>();
        for (Json varj : input["variations"].arrayValue())
            vars.push_back(std::shared_ptr<var_t>(var_t::parseVariation(varj)));
        Json colorj;
        if (color_dims && input.valueAt("color",colorj))
        {
            JsonArray colorja = colorj.arrayValue();
            if (colorja.size() != color_dims)
                throw std::runtime_error("color length incorrect");
            color.resize(color_dims);
            for (size_t i = 0; i < color_dims; ++i)
            {
                color[i] = colorja[i].floatValue();
                if (color[i] < 0.0 || color[i] > 1.0)
                    throw std::runtime_error("color coordinate out of range");
            }
        }
        _optimize();
    }
    inline num_t getWeight() const
    {
        return weight;
    }
    inline size_t getID() const
    {
        return id;
    }
    inline const Affine<num_t,dims>& getPreAffine() const
    {
        return pre;
    }
    inline const Affine<num_t,dims>& getPostAffine() const
    {
        return post;
    }
    inline const std::vector<std::shared_ptr<var_t>>& getVariations() const
    {
        return vars;
    }
    // iterate a state for the rendering process
    inline point_t applyIteration(rng_t<num_t>& rng, const point_t& p) const
    {
        point_t t,v;
        // for 2d, faster to apply identity affine than to branch
        if (dims < 3 || has_pre)
            t = pre.apply_to(p);
        for (auto var : vars)
            v += var->getWeight() * var->calc(rng,t);
        if (dims < 3 || has_post)
            t = post.apply_to(v);
        return t;
    }
    inline const std::vector<num_t>& getColor() const
    {
        return color;
    }
    inline bool hasColor() const
    {
        return !color.empty();
    }
};

}
