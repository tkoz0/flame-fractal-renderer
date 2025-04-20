#pragma once

#include "render_iterator.hpp"

#include "../types/constants.hpp"
#include "../types/flame.hpp"
#include "../types/types.hpp"

#include <atomic>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

namespace tkoz::flame
{

// union type for the buffer
union buf_elem_t
{
    hist_t uintval;
    num_t floatval;
    std::atomic<hist_t> uintval_atomic;
    std::atomic<num_t> floatval_atomic;
    // need to define constructor for C++20
    buf_elem_t(){}
};
static_assert(sizeof(buf_elem_t) == sizeof(hist_t));

/*
Buffer Renderer
- renders a buffer on a rectangular region
- includes histogram
- includes color information if enabled and flame includes color
- kept simple to leave the complex visualization stuff to an image renderer
*/
template <size_t dims, bool enable_color = true>
class BufferRenderer
{
    static_assert(dims > 0);
    static_assert(sizeof(num_t) == sizeof(hist_t));

    private:

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

    // flame to be rendered
    flame_t flame;
    // the buffer, interpreted as an array of num_t and hist_t
    std::vector<buf_elem_t> buffer;
    // buffer length in cells (1 per rendering region)
    size_t buffer_cells;
    // size of a cell (1 for the histogram + 1 for each color dimension)
    size_t buffer_cell_size;

    // get pointer to first element of a buffer cell
    [[nodiscard]] inline buf_elem_t *_buf_cell(size_t i)
    {
        if (enable_color)
            return buffer.data() + (i*buffer_cell_size);
        else
            return buffer.data() + i;
    }

    // get const pointer to first element of a buffer cell
    [[nodiscard]] inline const buf_elem_t *_buf_cell(size_t i) const
    {
        if (enable_color)
            return buffer.data() + (i*buffer_cell_size);
        else
            return buffer.data() + i;
    }

    // rendering statistics
    struct
    {
        // mutex for controlling thread access
        std::mutex mutex;
        // samples iterated during rendering
        size_t s_iter;
        // samples plotted during rendering
        size_t s_plot;
        // xform distribution by id
        std::vector<size_t> xf_dist;
        // bad value xform ids
        std::vector<size_t> bv_xfs;
        // bad value points
        std::vector<point_t> bv_pts;
        // min/max coordinates reached during iteration
        std::array<num_pair_t,dims> pt_max;
        // lock the mutex
        inline void lock() { mutex.lock(); }
        // unlock the mutex
        inline void unlock() { mutex.unlock(); }
    }
    stats;

    // multipliers for index in a single dimension
    std::array<num_t,dims> mult_d;
    // multipliers for calculating indexes
    std::array<size_t,dims> mult_i;

    // shared constructor stuff for initializing things from the flame params
    void _init()
    {
        stats.s_iter = 0;
        stats.s_plot = 0;
        stats.xf_dist = std::vector<size_t>(flame.getXFormIDCount(),0);
        // buffer size setup
        buffer_cells = 1;
        buffer_cell_size = 1;
        if (enable_color)
            buffer_cell_size += flame.getColorDims();
        for (size_t i = 0; i < dims; ++i)
        {
            size_t size = flame.getSize()[i];
            const num_pair_t& bounds = flame.getBounds()[i];
            mult_d[i] = (num_t)(size) / (bounds.second - bounds.first);
            mult_d[i] *= scale_adjust_down_v<num_t>;
            mult_i[i] = buffer_cells;
            buffer_cells *= size;
            if (buffer_cells >= (1uLL << 48))
                throw std::runtime_error("BufferRenderer(): histogram too big");
            stats.pt_max[i].first = INFINITY;
            stats.pt_max[i].second = -INFINITY;
        }
        buffer = std::vector<buf_elem_t>(buffer_cells*buffer_cell_size);
        std::for_each(buffer.begin(),buffer.end(),
            [](buf_elem_t& b) { b.uintval = 0; });
    }

    /*
    render a batch into the buffer
    - multithreaded = handle safe access to shared resources
    - samples = number of samples to render
    - bv_limit = number of bad/diverging points before terminating render
    - rng = random number generator object to use
    returns false if terminating due to bad value limit
    */
    template <bool multithreaded>
    bool _render_batch(size_t samples, size_t bv_limit)
    {
        if (stats.bv_xfs.size() > bv_limit)
            return false;
        RenderIterator<dims,enable_color> iter(flame);
        // local stats variables to avoid locking mutex on every iteration
        num_t s_iter = 0;
        num_t s_plot = 0;
        size_t *xf_dist = new size_t[flame.getXFormIDCount()]();
        std::array<num_pair_t,dims> pt_max;
        if (multithreaded)
            stats.lock();
        pt_max = stats.pt_max;
        if (multithreaded)
            stats.unlock();
        // begin render
        for (; samples > 0; --samples)
        {
            size_t xf_id;
            iter.iterate(xf_id);
            ++s_iter;
            ++xf_dist[xf_id];
            const point_t& p = iter.getPoint();
            // bad value check
            if (iter.badValue()) [[unlikely]]
            {
                if (multithreaded)
                    stats.lock();
                stats.bv_xfs.push_back(xf_id);
                stats.bv_pts.push_back(p);
                if (multithreaded)
                    stats.unlock();
                if (stats.bv_xfs.size() > bv_limit) [[unlikely]]
                    break;
                iter.init();
            }
            // extreme coordinates
            for (size_t i = 0; i < dims; ++i)
            {
                if (p[i] < pt_max[i].first) [[unlikely]]
                    pt_max[i].first = p[i];
                if (p[i] > pt_max[i].second) [[unlikely]]
                    pt_max[i].second = p[i];
            }
            // plot point
            if (!iter.inBoundsFinal())
                continue;
            ++s_plot;
            const point_t& pf = iter.getPointFinal();
            const std::array<num_pair_t,dims>& bounds = flame.getBounds();
            // calculate buffer index, starting with dimension 0 index
            size_t bi = (pf[0] - bounds[0].first) * mult_d[0];
            for (size_t i = 1; i < dims; ++i)
            {
                // calculate dimension i index
                size_t di = (pf[i] - bounds[i].first) * mult_d[i];
                // add to buffer index
                bi += di * mult_i[i];
            }
            // histogram increment
            buf_elem_t *bptr = _buf_cell(bi);
            if (multithreaded)
                ++(bptr->uintval_atomic);
            else
                ++(bptr->uintval);
            // add color
            if (enable_color)
            {
                const num_t *cf = iter.getColorFinal();
                size_t r = iter.getColorDims();
                if (multithreaded)
                    // fetch_add for float/double requires C++20
                    for (size_t i = 0; i < r; ++i)
                        (++bptr)->floatval_atomic.fetch_add(cf[i],
                            std::memory_order_relaxed);
                else
                    for (size_t i = 0; i < r; ++i)
                        (++bptr)->floatval += cf[i];
            }
        }
        // update statistics
        if (multithreaded)
            stats.lock();
        stats.s_iter += s_iter;
        stats.s_plot += s_plot;
        for (size_t i = 0; i < flame.getXFormIDCount(); ++i)
            stats.xf_dist[i] += xf_dist[i];
        for (size_t i = 0; i < dims; ++i)
        {
            stats.pt_max[i].first =
                std::min(stats.pt_max[i].first,pt_max[i].first);
            stats.pt_max[i].second =
                std::max(stats.pt_max[i].second,pt_max[i].second);
        }
        if (multithreaded)
            stats.unlock();
        // cleanup
        delete[] xf_dist;
        return samples == 0;
    }

public:

    BufferRenderer(const flame_t& flame): flame(flame) { _init(); }

    BufferRenderer(const Json& json): flame(json) { _init(); }

    /*
    render samples into the buffer
    - num_samples = number of samples to render
    - num_threads = maximum number of threads to use (capped to 65535)
    - batch_size = number of samples per work unit
    - bv_limit = number of bad/diverging points before terminating render
    - cb_batch = callback with progress on batch completion
      - called with a number in [0,1] indicating rendering progress
    - cb_thread = callback on thread start
      - called with a thread reference and a 0-indexed thread number
    */
    bool render(size_t num_samples, size_t num_threads, size_t batch_size,
            size_t bv_limit, std::function<void()> cb_batch = nullptr,
            std::function<void(const std::thread&,size_t)> cb_thread = nullptr)
    {
        stats.bv_pts.clear();
        stats.bv_xfs.clear();
        if (num_samples == 0)
            return true;
        if (num_threads == 0)
            throw std::runtime_error(
                "BufferRenderer::render(): threads must be positive");
        if (num_threads > 65535)
            throw std::runtime_error(
                "BufferRenderer::render(): too many threads");
        if (batch_size < 256)
            throw std::runtime_error(
                "BufferRenderer::render(): batch size too small");
        // total number of batches to render
        size_t num_batches = (num_samples+batch_size-1) / batch_size;
        // do not use more threads than batches
        num_threads = std::min(num_threads,num_batches);
        // batch request lock
        std::mutex batch_lock;
        // total number of samples requested by threads to render
        size_t s_render = 0;
        // did all batches complete without reaching bad value limit
        bool success = true;
        // function for each thread to use
        auto batch_processor = [this,&batch_lock,&s_render,&success,
            &num_samples,&batch_size,&bv_limit,&cb_batch]()
        {
            rng::setSeed(); // seed rng for life of this thread
            for (;;) // loop requesting batches until done
            {
                batch_lock.lock();
                size_t render_size = std::min(batch_size,num_samples-s_render);
                s_render += render_size;
                batch_lock.unlock();
                // all batches completed
                if (render_size == 0) [[unlikely]]
                    break;
                // make sure batch completes successfully
                if (!_render_batch<true>(render_size,bv_limit)) [[unlikely]]
                {
                    success = false;
                    break;
                }
                if (cb_batch)
                {
                    batch_lock.lock();
                    cb_batch();
                    batch_lock.unlock();
                }
            }
        };
        // initialize threads, start all before first batch begins
        std::vector<std::thread> threads;
        batch_lock.lock();
        for (size_t i = 0; i < num_threads; ++i)
        {
            threads.push_back(std::thread(batch_processor));
            if (cb_thread)
                cb_thread(threads.back(),i);
        }
        batch_lock.unlock();
        // wait for all threads to finish
        std::for_each(threads.begin(),threads.end(),
            [](std::thread& t) { t.join(); });
        return success;
    }

    /*
    render samples into buffer deterministically
    uses a single thread and provided rng object
    - samples = number of samples to render
    - batch_size = number of samples per work unit (completed sequentially)
    - bv_limit = number of bad/diverging points before terminating render
    - cb_batch = callback with progress on batch completion
      - called with a number in [0,1] indicating rendering progress
    */
    bool renderSeeded(size_t num_samples, size_t batch_size, size_t bv_limit,
            std::function<void()> cb_batch = nullptr)
    {
        stats.bv_pts.clear();
        stats.bv_xfs.clear();
        if (num_samples == 0)
            return true;
        if (batch_size == 0)
            throw std::runtime_error(
                "BufferRenderer::render(): batch size must be positive");
        size_t s_render = 0;
        while (s_render < num_samples)
        {
            size_t render_size = std::min(batch_size,num_samples-s_render);
            if (!_render_batch<false>(render_size,bv_limit)) [[unlikely]]
            {
                return false;
            }
            s_render += render_size;
            if (cb_batch)
                cb_batch();
        }
        return true;
    }

    // add another buffer from pointer, must match the same buffer format
    template <typename T>
    void addBuffer(const T *buf)
    {
        static_assert(sizeof(T) == sizeof(buf_elem_t));
        const buf_elem_t *src = (const buf_elem_t*) buf;
        buf_elem_t *dest = buffer.data();
        for (size_t i = 0; i < buffer_cells; ++i)
        {
            // add counter
            (dest++)->uintval += (src++)->uintval;
            // add color
            if (enable_color)
            {
                for (size_t j = 0; j < flame.getColorDims(); ++j)
                    (dest++)->floatval += (src++)->floatval;
            }
        }
    }

    // add another buffer from vector, must match the same buffer format
    template <typename T>
    void addBuffer(const std::vector<T>& buf)
    {
        static_assert(sizeof(T) == sizeof(buf_elem_t));
        if (buf.size() != buffer.size())
            throw std::runtime_error(
                "BufferRenderer::addBuffer(): sizes do not match");
        addBuffer(buf.data());
    }

    // add buffer from another renderer, must match the same buffer format
    void addBuffer(const BufferRenderer<dims,enable_color>& renderer)
    {
        if (buffer_cells != renderer.buffer_cells
         || buffer_cell_size != renderer.buffer_cell_size)
            throw std::runtime_error(
                "BufferRenderer::addBuffer(): formats do not match");
        addBuffer(renderer.buffer);
    }

    // add buffer from a file, must match the same buffer format
    // returns true if successful
    // the buffer may be partially overwritten if this fails
    // use alloc_tmp = true to avoid this by storing the input buffer first
    template <bool alloc_tmp = false>
    bool addBuffer(std::istream& is)
    {
        if (alloc_tmp)
        {
            std::vector<buf_elem_t> buf(buffer.size());
            is.read((char*)buf.data(),buf.size()*sizeof(hist_t));
            if (!is.good())
                return false;
            addBuffer(buf);
        }
        else
        {
            buf_elem_t *dest = buffer.data();
            const buf_elem_t *src;
            std::vector<buf_elem_t> buf(buffer_cell_size);
            for (size_t i = 0; i < buffer_cells; ++i)
            {
                is.read((char*)buf.data(),buf.size()*sizeof(hist_t));
                if (!is.good())
                    return false;
                src = buf.data();
                // add counter
                (dest++)->uintval += (src++)->uintval;
                // add color
                if (enable_color)
                {
                    for (size_t j = 0; j < flame.getColorDims(); ++j)
                        (dest++)->floatval += (src++)->floatval;
                }
            }
        }
        return true;
    }

    // overwrite the current buffer, return true if successful
    // the buffer may be partially overwritten if this fails
    // use alloc_tmp = true to avoid this by storing the input buffer first
    template <bool alloc_tmp = false>
    bool readBuffer(std::istream& is)
    {
        if (alloc_tmp)
        {
            std::vector<buf_elem_t> buf(buffer.size());
            is.read((char*)buf.data(),buf.size()*sizeof(hist_t));
            if (!is.good())
                return false;
            buffer = buf;
        }
        else
        {
            is.read((char*)buffer.data(),buffer.size()*sizeof(hist_t));
            return is.good();
        }
    }

    // write the buffer data to an output stream, return true if successful
    bool writeBuffer(std::ostream& os) const
    {
        os.write((char*)buffer.data(),buffer.size()*sizeof(hist_t));
        return os.good();
    }

    // the total number of samples rendered in this buffer
    [[nodiscard]] size_t histogramSum() const
    {
        size_t ret = 0;
        for (size_t i = 0; i < buffer_cells; ++i)
            ret += _buf_cell(i)->uintval;
        return ret;
    }

    // the minimum cell value in the buffer histogram
    [[nodiscard]] hist_t histogramMin() const
    {
        hist_t ret = -1;
        for (size_t i = 0; i < buffer_cells; ++i)
            ret = std::min(ret,_buf_cell(i)->uintval);
        return ret;
    }

    // the maximum cell value in the buffer histogram
    [[nodiscard]] hist_t histogramMax() const
    {
        hist_t ret = 0;
        for (size_t i = 0; i < buffer_cells; ++i)
            ret = std::max(ret,_buf_cell(i)->uintval);
        return ret;
    }

    [[nodiscard]] inline const flame_t& getFlame() const
    {
        return flame;
    }

    [[nodiscard]] inline const std::vector<buf_elem_t>& getBuffer() const
    {
        return buffer;
    }

    [[nodiscard]] inline const buf_elem_t *getBufferCell(size_t i) const
    {
        return _buf_cell(i);
    }

    [[nodiscard]] inline size_t getBufferNumCells() const
    {
        return buffer_cells;
    }

    [[nodiscard]] inline size_t getBufferCellSize() const
    {
        return buffer_cell_size;
    }

    [[nodiscard]] inline size_t getSamplesIterated() const
    {
        return stats.s_iter;
    }

    [[nodiscard]] inline size_t getSamplesPlotted() const
    {
        return stats.s_plot;
    }

    [[nodiscard]] inline
    const std::vector<size_t>& getXFormDistribution() const
    {
        return stats.xf_dist;
    }

    [[nodiscard]] inline const std::vector<size_t>& getBadValueXForms() const
    {
        return stats.bv_xfs;
    }

    [[nodiscard]] inline const std::vector<point_t>& getBadValuePoints() const
    {
        return stats.bv_pts;
    }

    [[nodiscard]] inline
    const std::array<num_pair_t,dims>& getPointExtremes() const
    {
        return stats.pt_max;
    }

    [[nodiscard]] inline const std::array<num_t,dims>& getDimMults() const
    {
        return mult_d;
    }

    [[nodiscard]] inline const std::array<size_t,dims>& getIndexMults() const
    {
        return mult_i;
    }

    [[nodiscard]] inline size_t getDims() const
    {
        return dims;
    }

    [[nodiscard]] inline size_t getColorDims() const
    {
        return enable_color ? flame.getColorDims() : 0;
    }
};

} // namnespace tkoz::flame
