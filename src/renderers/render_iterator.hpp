#pragma once

#include "../types/constants.hpp"
#include "../types/flame.hpp"
#include "../types/types.hpp"

namespace tkoz::flame
{

/*
Render Iterator
- iterates a sequence of points/colors

TODO
- implement xaos (from apophysis)
  - instead of a single weight parameter per xform
    use a matrix, xform selection probability depends on previous xform
- implement linked xform (from apophysis)
*/
template <size_t dims, bool enable_color = true>
class RenderIterator
{
private:

    static_assert(dims > 0);

    // point type
    typedef Point<num_t,dims> point_t;
    // base variation type
    typedef vars::Variation<dims> var_t;
    // xform type
    typedef XForm<dims> xform_t;
    // flame type
    typedef Flame<dims> flame_t;
    // number pair
    typedef std::pair<num_t,num_t> num_pair_t;

    // iterating point
    point_t p;
    // iterating point after final xform
    point_t pf;
    // iterating color
    num_t *c;
    // iterating color after final xform
    num_t *cf;
    // flame to be rendered
    const flame_t& flame;
    // color dimensions (or 0 if disabled)
    size_t r;

    // initialize the iterating state (point and possibly color)
    inline void _init()
    {
        p = rng::randPoint<dims>();
        for (size_t i = 0; i < settle_iters_v<num_t>; ++i)
            p = flame.getRandomXForm().applyIteration(p);
        if (enable_color)
            for (size_t i = 0; i < r; ++i)
                c[i] = rng::randNum();
    }

    // check point for bad value
    [[nodiscard]] inline static bool _bad_value(const point_t& x)
    {
        for (size_t i = 0; i < dims; ++i)
            if (bad_value(x[i])) [[unlikely]]
                return true;
        return false;
    }

    // check if point is within rectangle bounds
    [[nodiscard]] inline bool _in_bounds(const point_t& x) const
    {
        const std::array<num_pair_t,dims>& bounds = flame.getBounds();
        for (size_t i = 0; i < dims; ++i)
            if (x[i] < bounds[i].first || x[i] > bounds[i].second) [[unlikely]]
                return false;
        return true;
    }

public:

    RenderIterator(const flame_t& flame): flame(flame)
    {
        r = enable_color ? flame.getColorDims() : 0;
        if (enable_color && r > 0)
        {
            c = new num_t[r]();
            cf = new num_t[r]();
        }
        else
            c = cf = nullptr;
        _init();
    }

    ~RenderIterator()
    {
        if (c)
            delete[] c;
        if (cf)
            delete[] cf;
    }

    // perform an iteration to the next point (and possibly color)
    // store the selected xform id in xf_id
    inline void iterate(size_t& xf_id)
    {
        const xform_t& xf = flame.getRandomXForm();
        num_t s = xf.getColorSpeed();
        xf_id = xf.getID();
        p = xf.applyIteration(p);
        if (enable_color && xf.hasColor())
        {
            for (size_t i = 0; i < r; ++i)
                c[i] = (1.0-s)*c[i] + s*xf.getColor()[i];
        }
        if (flame.hasFinalXForm())
        {
            const xform_t& xff = flame.getFinalXForm();
            s = xff.getColorSpeed();
            pf = xff.applyIteration(p);
            if (enable_color)
            {
                if (xff.hasColor())
                    for (size_t i = 0; i < r; ++i)
                        cf[i] = (1.0-s)*c[i] + s*xff.getColor()[i];
                else
                    for (size_t i = 0; i < r; ++i)
                        cf[i] = c[i];
            }
        }
        else
        {
            pf = p;
            if (enable_color)
                for (size_t i = 0; i < r; ++i)
                    cf[i] = c[i];
        }
    }

    // perform an iteration to the next point (and possibly color)
    inline void iterate()
    {
        size_t xf_id;
        iterate(xf_id);
    }

    // initialize point again for starting
    inline void init()
    {
        _init();
    }

    // check if point has reached a bad value
    [[nodiscard]] inline bool badValue() const
    {
        return _bad_value(p);
    }

    // check if final point has reached a bad value
    [[nodiscard]] inline bool badValueFinal() const
    {
        return _bad_value(pf);
    }

    // check if point is in bounds
    [[nodiscard]] inline bool inBounds() const
    {
        return _in_bounds(p);
    }

    // check if final point is in bounds
    [[nodiscard]] inline bool inBoundsFinal() const
    {
        return _in_bounds(pf);
    }

    [[nodiscard]] inline const point_t& getPoint() const
    {
        return p;
    }

    [[nodiscard]] inline const point_t& getPointFinal() const
    {
        return pf;
    }

    // should not use this if color is disabled
    [[nodiscard]] inline const num_t *getColor() const
    {
        return c;
    }

    // should not use this if color is disabled
    [[nodiscard]] inline const num_t *getColorFinal() const
    {
        return cf;
    }

    [[nodiscard]] inline const flame_t& getFlame() const
    {
        return flame;
    }

    [[nodiscard]] inline size_t getColorDims() const
    {
        return r;
    }
};

} // namespace tkoz::flame
