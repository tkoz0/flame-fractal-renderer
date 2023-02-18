/*
Representation of an xform
*/

#pragma once

#include <cstdlib>
#include <vector>

// forward declaration
namespace tkoz::flame
{
template <typename num_t, size_t dims> class XForm;
}

#include "../settings.hpp"
#include "point.hpp"
#include "affine.hpp"
#include "iter_state.hpp"
#include "../variations.hpp"

namespace tkoz::flame
{

// xform (including final xform)
template <typename num_t, size_t dims>
class XForm
{
private:
    num_t weight; // xform probability weight, not applicable for final xform
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
            post = Affine<num_t,dims>();
        for (Json varj : input["variations"].arrayValue())
            vars.push_back(std::shared_ptr<var_t>(var_t::parseVariation(varj)));
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
    inline const std::vector<std::shared_ptr<var_t>>& getVariations() const
    {
        return vars;
    }
    // iterate a state for the rendering process
    inline void applyIteration(IterState<num_t,dims>& state) const
    {
        // for 2d, faster to apply identity affine than to branch
        if (dims < 3 || has_pre)
            state.t = pre.apply_to(state.p);
        state.v = Point<num_t,dims>();
        for (auto var : vars)
            state.v += var->getWeight() * var->calc(*state.rng,state.t);
        if (dims < 3 || has_post)
            state.p = post.apply_to(state.v);
    }
};

#if INSTANTIATE_TEMPLATES
extern template class XForm<float,1>;
extern template class XForm<float,2>;
extern template class XForm<float,3>;
extern template class XForm<float,4>;
extern template class XForm<float,5>;

extern template class XForm<double,1>;
extern template class XForm<double,2>;
extern template class XForm<double,3>;
extern template class XForm<double,4>;
extern template class XForm<double,5>;
#endif

}
