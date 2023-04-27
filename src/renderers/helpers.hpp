#pragma once

#include "histogram_renderer.hpp"

namespace tkoz::flame
{

template <typename num_t, size_t dims, typename hist_t, bool use_cache>
struct _ihr_helper
{
    static_assert(dims > 1);
    static HistogramRendererInterface<num_t,hist_t,use_cache>*
    helper(const Json& j)
    {
        if (j["dimensions"].intValue() == dims)
            return new HistogramRenderer<num_t,dims,hist_t,use_cache>(j);
        return _ihr_helper<num_t,dims-1,hist_t,use_cache>::helper(j);
    }
};

template <typename num_t, typename hist_t, bool use_cache>
struct _ihr_helper<num_t,1,hist_t,use_cache>
{
    static HistogramRendererInterface<num_t,hist_t,use_cache>*
    helper(const Json& j)
    {
        if (j["dimensions"].intValue() == 1)
            return new HistogramRenderer<num_t,1,hist_t,use_cache>(j);
        return nullptr;
    }
};

// instantiate an instance of a HistogramRenderer based on number of dimensions
// returns a pointer to newly allocated memory or nullptr if failure
template <typename num_t, typename hist_t, bool use_cache = false>
HistogramRendererInterface<num_t,hist_t,use_cache>*
instantiateHistogramRenderer(const Json& j)
{
    // the 5 means allow 1,2,3,4,5 dimensions only
    return _ihr_helper<num_t,5,hist_t,use_cache>::helper(j);
}

}
