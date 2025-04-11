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
template <size_t dims> class XForm;
}

#include "point.hpp"
#include "affine.hpp"
#include "../variations/variations.hpp"

namespace tkoz::flame
{

// xform (including final xform)
template <size_t dims>
class XForm
{
private:

    typedef Point<num_t,dims> point_t;
    typedef vars::Variation<dims> var_t;
    // xform id unique to the flame (using order of input file)
    size_t id;
    // xform probability weight, not applicable for final xform
    num_t weight;
    // variations list (polymorphic)
    std::vector<std::shared_ptr<var_t>> vars;
    // pre affine
    Affine<num_t,dims> pre;
    bool has_pre;
    // post affine
    Affine<num_t,dims> post;
    bool has_post;
    // color vector, empty means no color
    std::vector<num_t> color;
    // color speed
    num_t color_speed;

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
    // input = json data to parse
    // id = 0-indexed xform id
    // is_final = is this a final xform
    // color_dims = dimension of color space
    XForm(const Json& input, size_t id, bool is_final, size_t color_dims,
            num_t default_color_speed)
    {
        this->id = id;
        if (!is_final)
        {
            try
            {
                weight = input["weight"].floatValue();
            }
            catch (std::exception& e)
            {
                throw JsonError("XForm(): cannot parse weight: "
                    + std::string(e.what()));
            }
        }
        else
            weight = 1.0; // unused
        if (weight < 0.0)
            throw JsonError("XForm(): weight is negative");
        Json affine;
        has_pre = input.valueAt("pre_affine",affine);
        if (has_pre)
        {
            try
            {
                pre = Affine<num_t,dims>(affine);
            }
            catch (std::exception& e)
            {
                throw JsonError("XForm(): cannot parse pre_affine"
                    + std::string(e.what()));
            }
        }
        else
            pre = Affine<num_t,dims>();
        has_post = input.valueAt("post_affine",affine);
        if (has_post)
        {
            try
            {
                post = Affine<num_t,dims>(affine);
            }
            catch (std::exception& e)
            {
                throw JsonError("XForm(): cannot parse post_affine"
                    + std::string(e.what()));
            }
        }
        else
            post = Affine<num_t,dims>();
        try
        {
            for (Json varj : input["variations"].arrayValue())
                vars.push_back(
                    std::shared_ptr<var_t>(var_t::parseVariation(varj)));
        }
        catch (std::exception& e)
        {
            throw JsonError("XForm(): cannot parse variations: "
                + std::string(e.what()));
        }
        Json colorj;
        if (color_dims && input.valueAt("color",colorj))
        {
            JsonArray colorja = colorj.arrayValue();
            try
            {
                if (colorja.size() != color_dims)
                    throw JsonError("color length incorrect");
                color.resize(color_dims);
                for (size_t i = 0; i < color_dims; ++i)
                {
                    color[i] = colorja[i].floatValue();
                    if (color[i] < 0.0 || color[i] > 1.0)
                        throw JsonError("color coordinate out of range");
                }
            }
            catch (std::exception& e)
            {
                throw JsonError("XForm(): cannot parse color: "
                    + std::string(e.what()));
            }
        }
        if (color_dims && input.valueAt("color_speed",colorj))
        {
            try
            {
                color_speed = colorj.floatValue();
                if (color_speed < 0.0 || color_speed > 1.0)
                    throw JsonError("color speed out of range");
            }
            catch (std::exception& e)
            {
                throw JsonError("XForm(): color speed issue: "
                    + std::string(e.what()));
            }
        }
        else
            color_speed = default_color_speed;
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

    inline bool hasPreAffine() const
    {
        return has_pre;
    }

    inline bool hasPostAffine() const
    {
        return has_post;
    }

    // iterate a state for the rendering process
    inline point_t applyIteration(rng_t& rng, const point_t& p) const
    {
        point_t t,v;
        // for 2d, faster to apply identity affine than to branch
        if (dims < 3 || has_pre)
            t = pre.apply_to(p);
        else
            t = p;
        for (auto var : vars)
            v += var->getWeight() * var->calc(rng,t);
        if (dims < 3 || has_post)
            t = post.apply_to(v);
        else
            t = v;
        return t;
    }

    inline const std::vector<num_t>& getColor() const
    {
        return color;
    }

    inline bool hasColor() const
    {
        return color.size() != 0;
    }

    inline num_t getColorSpeed() const
    {
        return color_speed;
    }
};

} // namespace tkoz::flame
