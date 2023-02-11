/*
Representation of an xform
*/

#pragma once

#include <cstdlib>
#include <vector>

#include "affine.hpp"
#include "../variations.hpp"

namespace tkoz
{
namespace flame
{

// forward declaration
// need this in the xform constructor
//template <typename num_t, size_t dims> struct Variations;

// xform (including final xform)
template <typename num_t, size_t dims>
class XForm
{
private:
    num_t weight; // xform probability weight, not applicable for final xform
//    std::vector<XFormVar<num_t,dims>> vars; // variations
//    std::vector<num_t> varp; // variation parameters
    typedef vars::Variation<num_t,dims> var_t;
    std::vector<std::shared_ptr<var_t>> vars;
    Affine<num_t,dims> pre; // pre affine transformation
    Affine<num_t,dims> post; // post affine transformation
    bool has_pre,has_post;
public:
    XForm(){}
    // construct from JSON data
    // throws error if something goes wrong
    XForm(const Json& input, bool is_final = false)
    {
        if (!is_final)
            weight = input["weight"].floatValue();
        else
            weight = 1.0; // unused
        if (weight <= 0.0)
            throw std::runtime_error("weights must be positive");
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
            post = Affine<num_t,2>();
        for (Json varj : input["variations"].arrayValue())
        {
            vars.push_back(std::shared_ptr<var_t>(var_t::parseVariation(varj)));
            //XFormVar<num_t,dims> var;
            //std::string name = varj["name"].stringValue();
            //const VarInfo<num_t,dims>& varinfo =
            //    Variations<num_t,dims>::get(name);
            //var.func = varinfo.getFPtr();
            //var.index = varp.size();
            //vars.push_back(var);
            //varinfo.getPPtr()(varj,varp);
        }
    }
    // optimize xform
    void optimize()
    {
    }
    inline num_t getWeight() const
    {
        return weight;
    }
    inline const Affine<num_t,dims>& getPreAffine() const
    {
        return pre;
    }
    inline const Affine<num_t,dims>& getPostAffine() const
    {
        return post;
    }
    inline const std::vector<var_t>& getVariations() const
    {
        return vars;
    }
    //inline const std::vector<num_t>& getVariationParams() const
    //{
    //    return varp;
    //}
    // iterate a state for the rendering process
    inline void applyIteration(IterState<num_t,dims>& state) const
    {
        // for 2d, faster to apply identity affine than to branch
        if (dims < 3 || has_pre)
            state.t = pre.apply_to(state.p);
        state.v = Point<num_t,dims>();
        for (auto var : vars)
            //v.func(state,varp.data()+v.index);
            var->calc(state);
        if (dims < 3 || has_post)
            state.p = post.apply_to(state.v);
    }
};

}
}
