#pragma once

#include <array>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "../types/types.hpp"
#include "../types/flame.hpp"
#include "../utils/hardware.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

namespace tkoz::flame
{

// type generic interface to allow runtime template selection
// without duplicating a lot of code
template <typename num_t, typename hist_t, bool use_cache = false>
class HistogramRendererInterface
{
public:
    virtual void render(size_t samples, size_t threads, size_t batch_size,
        size_t bv_limit, std::function<void(float)> cb_batch = nullptr,
        std::function<void(const std::thread&,size_t)> cb_thread = nullptr) = 0;
    virtual inline size_t getDimensions() const = 0;
    virtual bool writeHistogram(std::ostream& os) const = 0;
    virtual inline size_t getSamplesIterated() const = 0;
    virtual inline size_t getSamplesPlotted() const = 0;
    virtual inline size_t getBadValueCount() const = 0;
    virtual inline size_t getHistogramSize() const = 0;
    virtual inline size_t getHistogramSizeBytes() const = 0;
    virtual void addHistogram(const hist_t *buf) = 0;
};

/*
Histogram Renderer
- renders a histogram on a rectangular region
- probability distribution only, no color information

TODO
- consider different data structures
- possibly u16 for histogram but store overflowed numbers in a separate map
- output for counting histogram entries in different ranges
  - how many need 8 bit, 16 bit, ... integer size
*/
template <typename num_t, size_t dims, typename hist_t, bool use_cache = false>
class HistogramRenderer:
    public HistogramRendererInterface<num_t,hist_t,use_cache>
{
    static_assert(!use_cache || sizeof(hist_t) > sizeof(u8));
private:
    typedef Point<num_t,dims> point_t;
    typedef XForm<num_t,dims> xform_t;
    typedef std::pair<num_t,num_t> num_pair_t;
    // flame to render
    Flame<num_t,dims> flame;
    // multipliers used by renderer to calculate indexes
    std::array<num_t,dims> xmults;
    std::array<size_t,dims> dsizes;
    // histogram arrays
    hist_t *histogram;
    u8 *histcache;
    size_t histsize;
    // === render statistics ===
    std::mutex lock_stats;
    // count samples inside the rendering bounds
    size_t samples_iterated, samples_plotted;
    // distribution of xforms selected
    std::vector<hist_t> xf_freq;
    // store xform index and point when reaching a bad value
    std::vector<size_t> bv_xforms;
    std::vector<point_t> bv_points;
    // min/max of coordinates reached during render
    std::array<num_pair_t,dims> point_extremes;
    // no default constructor, require constructor with flame
    HistogramRenderer(){}
    // initialize point with the initial convergence iterations
    inline void _init_point(rng_t<num_t>& rng, point_t& p) const
    {
        p = rng.template randPoint<dims>();
        for (size_t i = 0; i < settle_iters<num_t>::value; ++i)
            p = flame.getRandomXForm(rng).applyIteration(rng,p);
    }
    // get value at histogram index, handling cache
    inline hist_t _hist_index(size_t i) const
    {
        hist_t ret = histogram[i];
        if (use_cache)
            ret += histcache[i];
        return ret;
    }
    /*
    render samples into the histogram
    - samples = number of samples to render
    - bv_limit = stop after this many bad values
    - rng = random number generator object
    - possible defaults:
      - samples = 100*getHistogramSize()
      - bv_limit = 10
      - rng = random initialized by time
    */
    void _render_thread(size_t samples, size_t bv_limit, rng_t<num_t>& rng)
    {
        if (bv_xforms.size() >= bv_limit)
            return;
        point_t p,pf;
        // local statistics to avoid synchronization issues
        size_t samples_iterated_local = 0;
        size_t samples_plotted_local = 0;
        hist_t *xf_freq_local = new hist_t[flame.getXFormIDCount()]();
        // rendering loop
        _init_point(rng,p);
        for (size_t s = 0; s < samples; ++s)
        {
            render_loop_1: // start of the main render loop
            const xform_t& xf = flame.getRandomXForm(rng);
            p = xf.applyIteration(rng,p);
            ++samples_iterated_local;
            ++xf_freq_local[xf.getID()];
            // check for bad values
            for (size_t i = 0; i < dims; ++i)
                if (unlikely(bad_value(p[i])))
                {
                    lock_stats.lock();
                    bv_xforms.push_back(xf.getID());
                    bv_points.push_back(p);
                    lock_stats.unlock();
                    // terminate render if too many bad values
                    if (bv_xforms.size() > bv_limit)
                        goto render_loop_2;
                    _init_point(rng,p);
                    ++s;
                    // using goto is normally not good practice but it is sort
                    // of needed here with the extra loop for dimensions
                    goto render_loop_1;
                }
            // update extreme coordinates
            for (size_t i = 0; i < dims; ++i)
            {
                if (unlikely(p[i] < point_extremes[i].first))
                    point_extremes[i].first = p[i];
                if (unlikely(p[i] > point_extremes[i].second))
                    point_extremes[i].second = p[i];
            }
            // final xform
            if (flame.hasFinalXForm())
                pf = flame.getFinalXForm().applyIteration(rng,p);
            else
                pf = p;
            // check bounds before plotting
            for (size_t i = 0; i < dims; ++i)
                if (pf[i] < flame.getBounds()[i].first
                        || pf[i] > flame.getBounds()[i].second)
                {
                    ++s;
                    goto render_loop_1; // return to start of render loop
                }
            // calculate histogram index
            const auto& bounds = flame.getBounds();
            size_t index = (pf[0]-bounds[0].first)*xmults[0];
            for (size_t i = 1; i < dims; ++i)
                index += dsizes[i]*(size_t)((pf[i]-bounds[i].first)*xmults[i]);
#if 1
            // histogram increment with atomic operations
            if (use_cache)
            {
                if (unlikely(!__atomic_add_fetch(histcache+index,1,
                        __ATOMIC_RELAXED)))
                    __atomic_fetch_add(histogram+index,256,__ATOMIC_RELAXED);
            }
            else
                __atomic_fetch_add(histogram+index,1,__ATOMIC_RELAXED);
#else
            // histogram increment without atomic operations
            if (use_cache)
            {
                if (unlikely(++histcache[index] == 0))
                    histogram[index] += 256;
            }
            else
                ++histogram[index];
#endif
            ++samples_plotted_local;
        }
        render_loop_2: // end of the render loop
        // update statistics
        lock_stats.lock();
        samples_iterated += samples_iterated_local;
        samples_plotted += samples_plotted_local;
        for (size_t i = 0; i < flame.getXForms().size(); ++i)
            xf_freq[i] += xf_freq_local[i];
        lock_stats.unlock();
        // cleanup
        delete[] xf_freq_local;
    }
public:
    HistogramRenderer(const Flame<num_t,dims>& flame): flame(flame)
    {
        samples_iterated = 0;
        samples_plotted = 0;
        histsize = 1;
        dsizes[0] = 1;
        for (size_t i = 0; i < dims; ++i)
        {
            size_t size = flame.getSize()[i];
            const num_pair_t& bounds = flame.getBounds()[i];
            xmults[i] = (num_t)size / (bounds.second - bounds.first);
            // correction for ensuring indexing is in bounds
            xmults[i] *= scale_adjust<num_t>::value;
            histsize *= size;
            if (histsize >= (1uLL << 48))
                throw std::runtime_error("histogram too big");
            point_extremes[i].first = INFINITY;
            point_extremes[i].second = -INFINITY;
            if (i < dims-1)
                dsizes[i+1] = dsizes[i] * size;
        }
        histogram = new hist_t[histsize]();
        if (use_cache)
            histcache = new u8[histsize]();
        xf_freq = std::vector<hist_t>(flame.getXFormIDCount(),0);
    }
    ~HistogramRenderer()
    {
        delete[] histogram;
        if (use_cache)
            delete[] histcache;
    }
    /*
    render samples into the histogram using multiple threads
    - samples = number of samples to render
    - threads = number of threads to use
    - batch_size = number of samples per work unit
    - bv_limit = stop after this many bad values
    - cb_batch = callback with progress every time a batch complete
    - cb_thread = callback each time a thread starts
    - possible defaults:
      - samples = 1000*getHistogramSize()
      - threads = number_of_threads()
      - batch_size = 1 << 16
      - bv_limit = 10
    */
    void render(size_t samples, size_t threads, size_t batch_size,
        size_t bv_limit, std::function<void(float)> cb_batch = nullptr,
        std::function<void(const std::thread&,size_t)> cb_thread = nullptr)
    {
        if (threads == 0)
            throw std::runtime_error("threads must be positive");
        if (batch_size == 0)
            throw std::runtime_error("batch size must be positive");
        size_t batches = samples/batch_size + (samples % batch_size > 0);
        threads = std::min(threads,batches);
        // work unit request state
        std::mutex lock_wu;
        size_t samples_left = samples;
        auto thread_func = [this,&samples,&threads,&samples_left,&lock_wu,
                &batch_size,&bv_limit,&cb_batch,&cb_thread]()
        {
            rng_t<num_t> rng;
            for (;;) // request work units until none left
            {
                lock_wu.lock();
                size_t batch_samples = std::min(samples_left,batch_size);
                samples_left -= batch_samples;
                lock_wu.unlock();
                if (batch_samples == 0)
                    break;
                _render_thread(batch_samples,bv_limit,rng);
                if (cb_batch)
                {
                    lock_wu.lock();
                    cb_batch((float)(samples-samples_left)/samples);
                    lock_wu.unlock();
                }
            }
        };
        std::vector<std::thread> tl;
        lock_wu.lock(); // ensure all threads start first
        for (size_t i = 0; i < threads; ++i)
        {
            tl.push_back(std::thread(thread_func));
            if (cb_thread)
                cb_thread(tl.back(),i);
        }
        lock_wu.unlock();
        std::for_each(tl.begin(),tl.end(),[](std::thread& t){ t.join(); });
    }
    /*
    render as grayscale image (only supported for 2d buffer currently)
    if buf is null, allocates one, otherwise use existing buffer
    returns the buffer storing the resulting image
    the scaling function must map onto nonnegative numbers
    the scaling function should be nondecreasing
    */
    template <typename pix_t>
    pix_t *renderImageBuffer(std::function<num_t(hist_t)> scale,
            pix_t *buf = nullptr)
    {
        if (dims != 2)
            throw std::runtime_error("not implemented");
        if (!buf)
            buf = new pix_t[histsize];
        size_t nx = flame.getSize()[0];
        size_t ny = flame.getSize()[1];
        // preprocess histogram
        hist_t min = -1, max = 0;
        num_t smin = INFINITY, smax = -INFINITY;
        for (size_t i = 0; i < histsize; ++i)
        {
            hist_t h = _hist_index(i);
            num_t sh = scale(h);
            min = std::min(min,h);
            max = std::max(max,h);
            smin = std::min(smin,sh);
            smax = std::max(smax,sh);
            // use this opportunity to add cache into the main histogram
            if (use_cache)
            {
                histogram[i] += histcache[i];
                histcache[i] = 0;
            }
        }
        num_t mult;
        if (smax <= 0.0) // render as all black image in this case
            mult = 0.0;
        else
            mult = pix_scale<pix_t,num_t>::value / smax;
        // fill buffer with scaled values
        pix_t *p = buf;
        for (size_t y = ny; y--;)
            for (size_t x = 0; x < nx; ++x)
            {
                hist_t h = histogram[x + nx*y];
                *(p++) = (pix_t)(scale(h)*mult);
            }
        return buf;
    }
    inline const Flame<num_t,dims>& getFlame() const
    {
        return flame;
    }
    inline size_t getSamplesIterated() const
    {
        return samples_iterated;
    }
    inline size_t getSamplesPlotted() const
    {
        return samples_plotted;
    }
    inline const std::vector<hist_t>& getXFormFrequency() const
    {
        return xf_freq;
    }
    inline const std::vector<size_t>& getBadValueXForms() const
    {
        return bv_xforms;
    }
    inline const std::vector<point_t>& getBadValuePoints() const
    {
        return bv_points;
    }
    inline size_t getBadValueCount() const
    {
        return bv_xforms.size();
    }
    inline const std::array<num_pair_t,dims>& getPointExtremes() const
    {
        return point_extremes;
    }
    inline size_t getHistogramSize() const
    {
        return histsize;
    }
    inline size_t getHistogramSizeBytes() const
    {
        return histsize * sizeof(hist_t);
    }
    // add another buffer array to the internal buffer
    // length must be equal to getHistogramSize()
    void addHistogram(const hist_t *buf)
    {
        for (size_t i = 0; i < histsize; ++i)
            histogram[i] += buf[i];
    }
    // number of samples rendered in the histogram
    size_t histogramSum() const
    {
        size_t ret = 0;
        for (size_t i = 0; i < histsize; ++i)
            ret += _hist_index(i);
        return ret;
    }
    // minimum value in the histogram
    hist_t histogramMin() const
    {
        hist_t ret = -1;
        for (size_t i = 0; i < histsize; ++i)
            ret = std::min(ret,_hist_index(i));
        return ret;
    }
    // maximum value in the histogram
    hist_t histogramMax() const
    {
        hist_t ret = 0;
        for (size_t i = 0; i < histsize; ++i)
            ret = std::max(ret,_hist_index(i));
        return ret;
    }
    // simultaneously compute sum,min,max
    void histogramStats(size_t& sum, hist_t& min, hist_t& max) const
    {
        sum = 0;
        min = -1;
        max = 0;
        for (size_t i = 0; i < histsize; ++i)
        {
            hist_t h = _hist_index(i);
            sum += h;
            min = std::min(min,h);
            max = std::max(max,h);
        }
    }
    // write buffer, returns whether it is successful
    bool writeHistogram(std::ostream& os) const
    {
        // add cache to main histogram
        if (use_cache)
            for (size_t i = 0; i < histsize; ++i)
            {
                histogram[i] += histcache[i];
                histcache[i] = 0;
            }
        os.write((char*)histogram,getHistogramSizeBytes());
        return os.good();
    }
    // get a pointer to the histogram
    // handle_cache specifies if cache needs to be combined first
    inline hist_t *getHistogram(bool handle_cache = true)
    {
        if (use_cache && handle_cache)
            for (size_t i = 0; i < histsize; ++i)
            {
                histogram[i] += histcache[i];
                histcache[i] = 0;
            }
        return histogram;
    }
    inline size_t getDimensions() const
    {
        return dims;
    }
};

}

#undef likely
#undef unlikely
