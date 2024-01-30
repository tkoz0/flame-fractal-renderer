#pragma once

#include "../types/constants.hpp"
#include "../types/flame.hpp"
#include "../types/types.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

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
    // some typedefs
    typedef Point<num_t,dims> point_t;
    typedef vars::Variation<dims> var_t;
    typedef XForm<dims> xform_t;
    typedef Flame<dims> flame_t;
    typedef std::pair<num_t,num_t> num_pair_t;
    // iterating point
    point_t p,pf;
    // iterating color
    num_t *c,*cf;
    // flame to be rendered
    const flame_t& flame;
    // color dimensions (or 0 if disabled)
    size_t r;
    // random number generator
    rng_t& rng;
    // initialize the iterating state (point and possibly color)
    inline void _init()
    {
        p = rng.randPoint<dims>();
        for (size_t i = 0; i < settle_iters<num_t>::value; ++i)
            p = flame.getRandomXForm(rng).applyIteration(rng,p);
        if (enable_color)
            for (size_t i = 0; i < r; ++i)
                c[i] = rng.randNum();
    }
    // check point for bad value
    inline static bool _bad_value(const point_t& x)
    {
        for (size_t i = 0; i < dims; ++i)
            if (unlikely(bad_value(x[i])))
                return true;
        return false;
    }
    // check if point is within rectangle bounds
    inline bool _in_bounds(const point_t& x) const
    {
        const std::array<num_pair_t,dims>& bounds = flame.getBounds();
        for (size_t i = 0; i < dims; ++i)
            if (x[i] < bounds[i].first || x[i] > bounds[i].second)
                return false;
        return true;
    }
public:
    RenderIterator(const flame_t& flame, rng_t& rng): flame(flame), rng(rng)
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
        const xform_t& xf = flame.getRandomXForm(rng);
        num_t s = xf.getColorSpeed();
        xf_id = xf.getID();
        p = xf.applyIteration(rng,p);
        if (enable_color && xf.hasColor())
        {
            for (size_t i = 0; i < r; ++i)
                c[i] = (1.0-s)*c[i] + s*xf.getColor()[i];
        }
        if (flame.hasFinalXForm())
        {
            const xform_t& xff = flame.getFinalXForm();
            s = xff.getColorSpeed();
            pf = xff.applyIteration(rng,p);
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
    inline bool badValue() const
    {
        return _bad_value(p);
    }
    // check if final point has reached a bad value
    inline bool badValueFinal() const
    {
        return _bad_value(pf);
    }
    // check if point is in bounds
    inline bool inBounds() const
    {
        return _in_bounds(p);
    }
    // check if final point is in bounds
    inline bool inBoundsFinal() const
    {
        return _in_bounds(pf);
    }
    inline const point_t& getPoint() const
    {
        return p;
    }
    inline const point_t& getPointFinal() const
    {
        return pf;
    }
    // should not use this if color is disabled
    inline const num_t *getColor() const
    {
        return c;
    }
    // should not use this if color is disabled
    inline const num_t *getColorFinal() const
    {
        return cf;
    }
    inline const flame_t& getFlame() const
    {
        return flame;
    }
    inline size_t getColorDims() const
    {
        return r;
    }
};

}

#undef likely
#undef unlikely
