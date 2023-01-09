#pragma once

#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

#include "types.hpp"
#include "variations.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

// macros for using SFINAE
#define FUNC_ENABLE_IFSAME(T1,T2,RET) template <typename dummy = T1> \
    typename std::enable_if<std::is_same<dummy,T2>::value,RET>::type
#define FUNC_ENABLE_IFSAME2(T1,T2,U1,U2,RET) template <typename dummy = T1> \
    typename std::enable_if<std::is_same<dummy,T2>::value \
        && std::is_same<U1,U2>::value,RET>::type

namespace tkoz
{
namespace flame
{

template <typename num_t, typename word_t, size_t rparam, size_t dims>
struct IterState<num_t,dims,Isaac<word_t,rparam>>
{
    Point<num_t,dims> p, t, v;
    Isaac<word_t,rparam>& rng;
    const XForm<num_t,2,Isaac<word_t,rparam>> *xf;
    num_t *cw;
    IterState(Isaac<word_t,rparam>& rng):
        p(),t(),v(),rng(rng),xf(nullptr) {}
    inline bool randBool() { return rng.next() & 1; }
    FUNC_ENABLE_IFSAME2(num_t,float,word_t,u32,num_t) inline randNum()
    { return (rng.next() >> 8) / (float)(1 << 24); }
    FUNC_ENABLE_IFSAME2(num_t,float,word_t,u64,num_t) inline randNum()
    { return (rng.next() >> 40) / (float)(1 << 24); }
    FUNC_ENABLE_IFSAME2(num_t,double,word_t,u32,num_t) inline randNum()
    {
        u32 hi = rng.next() >> 6;
        u32 lo = rng.next() >> 5;
        return (((u64)hi << 27) + lo) / (double)(1LL << 53);
    }
    FUNC_ENABLE_IFSAME2(num_t,double,word_t,u64,num_t) inline randNum()
    { return (rng.next() >> 11) / (double)(1LL << 53); }
    // slightly non-uniform distribution when mod is not a power of 2
    inline u32 randInt(u32 mod) { return (u32)(rng.next()) % mod; }
    // random point in the biunit square/cube/hypercube
    inline Point<num_t,dims> randPoint()
    {
        num_t x[dims];
        for (size_t i = 0; i < dims; ++i)
            x[i] = 2.0*randNum() - 1.0;
        return Point<num_t,dims>(x);
    }
    inline const Affine<num_t,2>& getPreAffine() const
    { return xf->getPreAffine(); }
    inline const Affine<num_t,2>& getPostAffine() const
    { return xf->getPostAffine(); }
    // use cumulative weights to select xform
    inline u32 randXFormIndex()
    {
        u32 ret = 0;
        num_t r = randNum();
        while (cw[ret] < r)
            ++ret;
        return ret;
    }
};

template <typename num_t, size_t dims, typename rand_t> struct XFormVar
{
    // function pointer
    std::function<void(IterState<num_t,dims,rand_t>&,const num_t*)> func;
    size_t index; // index of first variation parameter in varp (XForm class)
    // parameters are taken in order starting from index
    // the varp vector keeps the parameters compactly in memory
};

// xform (including final xform)
template <typename num_t, size_t dims, typename rand_t> class XForm
{
private:
    num_t weight; // xform probability weight, not applicable for final xform
    std::vector<XFormVar<num_t,dims,rand_t>> vars; // variations
    std::vector<num_t> varp; // variation parameters
    Affine<num_t,dims> pre; // pre affine transformation
    Affine<num_t,dims> post; // post affine transformation
    bool has_pre,has_post;
    u32 pc_flags; // precalculate flags
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
        //num_t A[6];
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
        pc_flags = 0;
        for (Json varj : input["variations"].arrayValue())
        {
            XFormVar<num_t,dims,rand_t> var;
            std::string name = varj["name"].stringValue();
            const VarInfo<num_t,dims,rand_t>& varinfo =
                Variations<num_t,dims,rand_t>::get(name);
            var.func = varinfo.func;
            var.index = varp.size();
            num_t weight = varj["weight"].floatValue();
            if (weight == 0.0)
                continue; // skip zero weight variations
            vars.push_back(var);
            if (varinfo.params) // store other parameters
                varinfo.params(*this,varj,weight,varp);
            else // store just the weight by default
                varp.push_back(weight);
            pc_flags |= varinfo.pc_flags;
        }
    }
    // optimize xform
    void optimize()
    {
        // zero weight variations are excluded in constructor
    }
    inline num_t getWeight() const { return weight; }
    inline const Affine<num_t,dims>& getPreAffine() const { return pre; }
    inline const Affine<num_t,dims>& getPostAffine() const { return post; }
    inline const std::vector<XFormVar<num_t,dims,rand_t>>& getVariations() const
    { return vars; }
    inline const std::vector<num_t>& getVariationParams() const
    { return varp; }
    // iterate a state for the rendering process
    inline void applyIteration(IterState<num_t,dims,rand_t>& state) const
    {
        state.xf = this;
        //if (has_pre)
            state.t = pre.apply_to(state.p);
        state.v = Point<num_t,dims>();
        for (XFormVar<num_t,dims,rand_t> v : vars)
            v.func(state,varp.data()+v.index);
        //if (has_post)
            state.p = post.apply_to(state.v);
    }
};

// flame fractal
template <typename num_t, size_t dims, typename rand_t> class Flame
{
private:
    size_t size_x,size_y; // dimensions
    num_t xmin,xmax,ymin,ymax; // rectangle bounds
    std::vector<XForm<num_t,dims,rand_t>> xforms;
    XForm<num_t,dims,rand_t> final_xform;
    bool has_final_xform;
    Flame(){}
public:
    // construct from JSON data
    // throws error if something goes wrong
    Flame(const Json& input)
    {
        // top level entries
        size_x = input["size_x"].intValue();
        size_y = input["size_y"].intValue();
        xmin = input["xmin"].floatValue();
        xmax = input["xmax"].floatValue();
        ymin = input["ymin"].floatValue();
        ymax = input["ymax"].floatValue();
        if (size_x == 0 || size_x > max_dim)
            throw std::runtime_error("size_x out of bounds");
        if (size_y == 0 || size_y > max_dim)
            throw std::runtime_error("size_y out of bounds");
        num_t M = max_rect<num_t>::value;
        if (xmin < -M || xmin > M)
            throw std::runtime_error("xmin out of bounds");
        if (xmax < -M || xmax > M)
            throw std::runtime_error("xmax out of bounds");
        if (ymin < -M || ymax > M)
            throw std::runtime_error("ymin out of bounds");
        if (ymax < -M || ymax > M)
            throw std::runtime_error("ymax out of bounds");
        if (xmin >= xmax)
            throw std::runtime_error("xmin >= xmax");
        if (ymin >= ymax)
            throw std::runtime_error("ymin >= ymax");
        // xforms loop
        for (Json xf : input["xforms"].arrayValue())
        {
            if (xf["weight"].floatValue() == 0.0) // ignore 0 weight xforms
                continue;
            xforms.push_back(XForm<num_t,dims,rand_t>(xf));
        }
        Json xf;
        has_final_xform = input.valueAt("final_xform",xf);
        if (has_final_xform)
            final_xform = XForm<num_t,dims,rand_t>(input["final_xform"],true);
    }
    void optimize()
    {
        // sort by decreasing weight
        std::sort(xforms.begin(),xforms.end(),
            [](XForm<num_t,dims,rand_t>& a, XForm<num_t,dims,rand_t>& b)
            { return a.getWeight() > b.getWeight(); });
        // optimize each xform
        std::for_each(xforms.begin(),xforms.end(),
            [](XForm<num_t,dims,rand_t>& xf) { xf.optimize(); });
    }
    inline size_t getSizeX() const { return size_x; }
    inline size_t getSizeY() const { return size_y; }
    inline num_t getXMin() const { return xmin; }
    inline num_t getXMax() const { return xmax; }
    inline num_t getYMin() const { return ymin; }
    inline num_t getYMax() const { return ymax; }
    inline const std::vector<XForm<num_t,dims,rand_t>>& getXForms() const
    { return xforms; }
    inline bool hasFinalXForm() const { return has_final_xform; }
    const XForm<num_t,dims,rand_t>& getFinalXForm() const { return final_xform; }
};

// render histogram only (count of samples in each pixel)
template <typename num_t, typename hist_t, typename rand_t>
class RendererBasic
{
private:
    Flame<num_t,2,rand_t> flame;
    hist_t *histogram;
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
public:
    // construct a renderer object from a flame, optionally an existing buffer
    // buf != null to use existing buffer, maybe loaded from a file
    RendererBasic(const Flame<num_t,2,rand_t>& flame, hist_t *buf = nullptr):
        flame(flame),samples_iterated(0),samples_plotted(0),
        xmin(INFINITY),ymin(INFINITY),
        xmax(-INFINITY),ymax(-INFINITY)
    {
        this->flame.optimize();
        hist_alloc = buf == nullptr;
        size_t X = flame.getSizeX();
        size_t Y = flame.getSizeY();
        if (!buf)
            histogram = new hist_t[X*Y]();
        else
            histogram = buf;
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
    }
    void renderBuffer(size_t samples, rand_t& rng, size_t bad_value_limit = 10)
    {
        if (bad_value_xforms.size() >= bad_value_limit)
            return;
        // state setup
        IterState<num_t,2,rand_t> state(rng);
        state.cw = cw;
        state.p = state.randPoint();
        const std::vector<XForm<num_t,2,rand_t>>& xfs = flame.getXForms();
        size_t flame_x = flame.getSizeX();
        size_t flame_y = flame.getSizeY();
        bool has_final_xform = flame.hasFinalXForm();
        // multipliers for calculating coordinates in histogram
        num_t xmul = (num_t)flame_x / (flame.getXMax() - flame.getXMin());
        num_t ymul = (num_t)flame_y / (flame.getYMax() - flame.getYMin());
        // correction to ensure indexing in bounds
        xmul *= scale_adjust<num_t>::value;
        ymul *= scale_adjust<num_t>::value;
        // get the point to converge to the attractor
        for (size_t s = 0; s < settle_iters<num_t>::value; ++s)
            xfs[state.randXFormIndex()].applyIteration(state);
        size_t samples_iterated_local = 0;
        size_t samples_plotted_local = 0;
        hist_t *xfdist_local = new hist_t[flame.getXForms().size()]();
        for (size_t s = 0; s < samples; ++s)
        {
            ++samples_iterated_local;
            u32 xf_i = state.randXFormIndex();
            const XForm<num_t,2,rand_t>& xf = xfs[xf_i];
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
                state.p = state.randPoint();
                for (size_t s = 0; s < settle_iters<num_t>::value; ++s)
                    xfs[state.randXFormIndex()].applyIteration(state);
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
            if (state.t.x() < flame.getXMin() || state.t.x() > flame.getXMax())
                continue;
            if (state.t.y() < flame.getYMin() || state.t.y() > flame.getYMax())
                continue;
            // increment in histogram
            size_t x = (state.t.x() - flame.getXMin()) * xmul;
            size_t y = (state.t.y() - flame.getYMin()) * ymul;
            //++histogram[flame_x*y + x];
            __atomic_fetch_add(histogram+(flame_x*y + x),1,__ATOMIC_RELAXED);
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
            rand_t rng; // unique random seeded rng for each thread
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
        for (size_t r = flame.getSizeY(); r--;)
            for (size_t c = 0; c < flame.getSizeX(); ++c)
            {
                hist_t buf_val = histogram[flame.getSizeX()*r + c];
                num_t scale_val = scale(buf_val);
                scale_min = std::min(scale_min,scale_val);
                scale_max = std::max(scale_max,scale_val);
            }
        pix_t *img_ptr = buf;
        num_t mult = pix_scale<pix_t,num_t>::value / scale_max;
        for (size_t r = flame.getSizeY(); r--;)
            for (size_t c = 0; c < flame.getSizeX(); ++c)
            {
                hist_t buf_val = histogram[flame.getSizeX()*r + c];
                *(img_ptr++) = (pix_t)(scale(buf_val)*mult);
            }
        return buf;
    }
    inline const Flame<num_t,2,rand_t>& getFlame() const { return flame; }
    inline const hist_t *getHistogram() const { return histogram; }
    inline hist_t *getHistogram() { return histogram; }
    size_t getHistogramSizeBytes() const
    { return sizeof(hist_t)*flame.getSizeX()*flame.getSizeY(); };
    size_t getHistogramSize() const
    { return flame.getSizeX()*flame.getSizeY(); }
    inline size_t getXFormsLength() const { return flame.getXForms().size(); }
    inline const hist_t *getXFormDistribution() const { return xfdist; }
    inline size_t getBadValueCount() const { return bad_value_xforms.size(); }
    inline size_t getSamplesPlotted() const { return samples_plotted; }
    inline size_t getSamplesIterated() const { return samples_iterated; }
    inline const std::vector<u32>& getBadValueXForms() const
    { return bad_value_xforms; }
    inline const std::vector<Point<num_t,2>>& getBadValuePoints() const
    { return bad_value_points; }
    inline num_t getXMin() const { return xmin; }
    inline num_t getXMax() const { return xmax; }
    inline num_t getYMin() const { return ymin; }
    inline num_t getYMax() const { return ymax; }
};

}
}

#undef FUNC_ENABLE_IF
#undef FUNC_ENABLE_IF2
#undef likely
#undef unlikely
