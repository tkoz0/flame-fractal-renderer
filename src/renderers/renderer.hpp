#pragma once

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

#include "../types/flame.hpp"
#include "../types/point.hpp"
#include "../types/types.hpp"
#include "../types/iter_state.hpp"
#include "../utils/flame.hpp"
#include "../utils/sfinae.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz::flame
{

// render histogram only (count of samples in each pixel)
template <typename num_t, typename hist_t>
class RendererBasic
{
    static_assert(sizeof(hist_t) > sizeof(u8));
private:
    Flame<num_t,2> flame;
    hist_t *histogram;
    // main write location for renderer threads (smaller fits in cache better)
    // on overflow, add to histogram
    u8 *bufcache;
    bool hist_alloc; // is histogram allocated by this instance
    num_t *cw; // cumulative weights for xform probability selection
    RendererBasic(){}
    // rendering statistics (current session)
    std::mutex mutex;
    size_t samples_iterated,samples_plotted;
    std::vector<u32> bad_value_xforms; // last xforms leading to bad value
    std::vector<Point<num_t,2>> bad_value_points; // points at bad value
    hist_t *xfdist; // xform selection (TODO maybe remove)
    // extreme coordinates during render, not handled atomically
    // only guaranteed to be correct with 1 thread
    num_t xmin,ymin,xmax,ymax;
    size_t X,Y;
public:
    // construct a renderer object from a flame, optionally an existing buffer
    // buf != null to use existing buffer, maybe loaded from a file
    RendererBasic(const Flame<num_t,2>& flame, hist_t *buf = nullptr):
        flame(flame),samples_iterated(0),samples_plotted(0),
        xmin(INFINITY),ymin(INFINITY),
        xmax(-INFINITY),ymax(-INFINITY),
        X(flame.getSize()[0]),Y(flame.getSize()[1])
    {
        this->flame.optimize();
        hist_alloc = buf == nullptr;
        if (!buf)
            histogram = new hist_t[X*Y]();
        else
            histogram = buf;
        bufcache = new u8[X*Y]();
        xfdist = new hist_t[flame.getXForms().size()]();
        // compute normalized xform weights
        std::vector<num_t> weights;
        num_t wsum = 0.0;
        for (auto xf : this->flame.getXForms())
        {
            weights.push_back(xf.getWeight());
            wsum += xf.getWeight();
        }
        for (size_t i = 0; i < weights.size(); ++i)
            weights[i] /= wsum;
        // cumulative weights for probability selection
        cw = new num_t[weights.size()];
        wsum = 0.0;
        for (size_t i = 0; i < weights.size(); ++i)
        {
            wsum += weights[i];
            cw[i] = wsum;
        }
        cw[weights.size()-1] = 1.0; // correct rounding error
    }
    ~RendererBasic()
    {
        // deallocate histogram if this instance allocated it
        if (hist_alloc)
            delete[] histogram;
        delete[] xfdist;
        delete[] cw;
        delete[] bufcache;
    }
    void renderBuffer(size_t samples, rng_t<num_t>& rng,
        size_t bad_value_limit = 10)
    {
        if (bad_value_xforms.size() >= bad_value_limit)
            return;
        // state setup
        IterState<num_t,2> state;
        state.rng = &rng;
        state.cw = cw;
        state.p = rng.template randPoint<2>();
        const std::vector<XForm<num_t,2>>& xfs = flame.getXForms();
        bool has_final_xform = flame.hasFinalXForm();
        // multipliers for calculating coordinates in histogram
        std::pair<num_t,num_t> xb = flame.getBounds()[0];
        std::pair<num_t,num_t> yb = flame.getBounds()[1];
        num_t xmul = (num_t)X / (xb.second - xb.first);
        num_t ymul = (num_t)Y / (yb.second - yb.first);
        // correction to ensure indexing in bounds
        xmul *= scale_adjust<num_t>::value;
        ymul *= scale_adjust<num_t>::value;
        // get the point to converge to the attractor
        for (size_t s = 0; s < settle_iters<num_t>::value; ++s)
            xfs[rng.randXFI(state.cw)].applyIteration(state);
        size_t samples_iterated_local = 0;
        size_t samples_plotted_local = 0;
        hist_t *xfdist_local = new hist_t[flame.getXForms().size()]();
        for (size_t s = 0; s < samples; ++s)
        {
            ++samples_iterated_local;
            size_t xf_i = rng.randXFI(state.cw);
            const XForm<num_t,2>& xf = xfs[xf_i];
            ++xfdist_local[xf_i];
            xf.applyIteration(state);
            if (unlikely(bad_value(state.p.x()) || bad_value(state.p.y())))
            {
                mutex.lock();
                bad_value_xforms.push_back(xf_i);
                bad_value_points.push_back(state.p);
                mutex.unlock();
                if (bad_value_xforms.size() >= bad_value_limit)
                    break;
                state.p = rng.template randPoint<2>();
                for (size_t s = 0; s < settle_iters<num_t>::value; ++s)
                    xfs[rng.randXFI(state.cw)].applyIteration(state);
                continue;
            }
            // update extreme coordinates
            if (unlikely(state.p.x() < xmin)) xmin = state.p.x();
            if (unlikely(state.p.x() > xmax)) xmax = state.p.x();
            if (unlikely(state.p.y() < ymin)) ymin = state.p.y();
            if (unlikely(state.p.y() > ymax)) ymax = state.p.y();
            if (has_final_xform) // update state.p to point to use
            {
                Point<num_t,2> tmp = state.p;
                flame.getFinalXForm().applyIteration(state);
                state.t = state.p;
                state.p = tmp;
            }
            else
                state.t = state.p;
            // skip plotting if out of bounds
            if (state.t.x() < xb.first || state.t.x() > xb.second)
                continue;
            if (state.t.y() < yb.first || state.t.y() > yb.second)
                continue;
            // increment in histogram
            size_t x = (state.t.x() - xb.first) * xmul;
            size_t y = (state.t.y() - yb.first) * ymul;
            size_t index = X*y + x;
            //++histogram[flame_x*y + x];
            //__atomic_fetch_add(histogram+(flame_x*y + x),1,__ATOMIC_RELAXED);
            // increment in cache, if overflow then add to histogram
            if (unlikely(!__atomic_add_fetch(bufcache+index,1,
                    __ATOMIC_RELAXED)))
                __atomic_fetch_add(histogram+index,1<<8,__ATOMIC_RELAXED);
            ++samples_plotted_local;
        }
        mutex.lock();
        samples_iterated += samples_iterated_local;
        samples_plotted += samples_plotted_local;
        for (size_t i = 0; i < flame.getXForms().size(); ++i)
            xfdist[i] += xfdist_local[i];
        mutex.unlock();
        delete[] xfdist_local;
    }
    void renderBufferParallel(size_t samples, size_t threads = 1,
            size_t batch_size = 1 << 16, size_t bad_value_limit = 10,
            std::function<void(float)> batch_callback = nullptr,
            std::function<void(std::thread&,size_t)> thread_callback = nullptr)
    {
        threads = std::min(threads,samples/batch_size+(samples%batch_size>0));
        std::mutex batch_mutex;
        size_t samples_progress = 0;
        size_t samples_total = samples;
        auto thread_function = [this,&samples,&batch_mutex,&samples_progress,
                &samples_total,&bad_value_limit,&batch_callback,&batch_size]()
        {
            rng_t<num_t> rng; // unique random seeded rng for each thread
            for (;;)
            {
                batch_mutex.lock(); // get next unit
                if (!samples)
                {
                    batch_mutex.unlock();
                    break;
                }
                size_t batch_samples = std::min(samples,batch_size);
                samples -= batch_samples;
                batch_mutex.unlock();
                renderBuffer(batch_samples,rng,bad_value_limit);
                batch_mutex.lock(); // completed unit
                samples_progress += batch_samples;
                if (batch_callback)
                    batch_callback((float)samples_progress/samples_total);
                batch_mutex.unlock();
            }
        };
        std::vector<std::thread> thread_list;
        batch_mutex.lock(); // start threads first
        for (size_t i = 0; i < threads; ++i)
        {
            thread_list.push_back(std::thread(thread_function));
            if (thread_callback)
                thread_callback(thread_list.back(),i);
        }
        batch_mutex.unlock();
        for (size_t i = 0; i < threads; ++i)
            thread_list[i].join();
        // add from cache to histogram
        for (size_t i = 0; i < getHistogramSize(); ++i)
        {
            histogram[i] += bufcache[i];
            bufcache[i] = 0;
        }
    }
    template <typename pix_t>
    pix_t *renderImage(std::function<num_t(hist_t)> scale, pix_t *buf = nullptr)
    {
        hist_t sample_min = -1;
        hist_t sample_max = 0;
        for (size_t i = 0; i < getHistogramSize(); ++i)
        {
            sample_min = std::min(sample_min,histogram[i]);
            sample_max = std::max(sample_max,histogram[i]);
        }
        if (!buf)
            buf = new pix_t[sizeof(pix_t)*getHistogramSize()];
        num_t scale_min = INFINITY;
        num_t scale_max = -INFINITY;
        for (size_t r = Y; r--;)
            for (size_t c = 0; c < X; ++c)
            {
                hist_t buf_val = histogram[X*r + c];
                num_t scale_val = scale(buf_val);
                scale_min = std::min(scale_min,scale_val);
                scale_max = std::max(scale_max,scale_val);
            }
        pix_t *img_ptr = buf;
        num_t mult = pix_scale<pix_t,num_t>::value / scale_max;
        for (size_t r = Y; r--;)
            for (size_t c = 0; c < X; ++c)
            {
                hist_t buf_val = histogram[X*r + c];
                *(img_ptr++) = (pix_t)(scale(buf_val)*mult);
            }
        return buf;
    }
    inline const Flame<num_t,2>& getFlame() const
    {
        return flame;
    }
    inline const hist_t *getHistogram() const
    {
        return histogram;
    }
    inline hist_t *getHistogram()
    {
        return histogram;
    }
    size_t getHistogramSizeBytes() const
    {
        return sizeof(hist_t)*X*Y;
    };
    size_t getHistogramSize() const
    {
        return X*Y;
    }
    inline size_t getXFormsLength() const
    {
        return flame.getXForms().size();
    }
    inline const hist_t *getXFormDistribution() const
    {
        return xfdist;
    }
    inline size_t getBadValueCount() const
    {
        return bad_value_xforms.size();
    }
    inline size_t getSamplesPlotted() const
    {
        return samples_plotted;
    }
    inline size_t getSamplesIterated() const
    {
        return samples_iterated;
    }
    inline const std::vector<u32>& getBadValueXForms() const
    {
        return bad_value_xforms;
    }
    inline const std::vector<Point<num_t,2>>& getBadValuePoints() const
    {
        return bad_value_points;
    }
    inline num_t getXMin() const
    {
        return xmin;
    }
    inline num_t getXMax() const
    {
        return xmax;
    }
    inline num_t getYMin() const
    {
        return ymin;
    }
    inline num_t getYMax() const
    {
        return ymax;
    }
};

}

#undef likely
#undef unlikely
