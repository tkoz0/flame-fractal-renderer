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
#define ENABLE_IF(COND,RET) template <typename RET2 = RET> \
    typename std::enable_if<(COND),RET2>::type
#define ENABLE_IFEQ(N1,N2,RET) template <typename RET2 = RET> \
    typename enable_if_eq<N1,N2,RET2>::type

namespace tkoz
{
namespace flame
{

template <typename num_t, typename word_t, size_t rparam, size_t dims>
struct IterState<num_t,dims,Isaac<word_t,rparam>>
{
    typedef Point<num_t,dims> point_t;
    point_t p, t, v;
    Isaac<word_t,rparam>& rng;
    const XForm<num_t,2,Isaac<word_t,rparam>> *xf;
    num_t *cw;
    IterState(Isaac<word_t,rparam>& rng):
        p(),t(),v(),rng(rng),xf(nullptr) {}
    inline bool randBool()
    {
        return rng.next() & 1;
    }
    FUNC_ENABLE_IFSAME2(num_t,float,word_t,u32,num_t) inline randNum()
    {
        return (rng.next() >> 8) / (float)(1 << 24);
    }
    FUNC_ENABLE_IFSAME2(num_t,float,word_t,u64,num_t) inline randNum()
    {
        return (rng.next() >> 40) / (float)(1 << 24);
    }
    FUNC_ENABLE_IFSAME2(num_t,double,word_t,u32,num_t) inline randNum()
    {
        u32 hi = rng.next() >> 6;
        u32 lo = rng.next() >> 5;
        return (((u64)hi << 27) + lo) / (double)(1LL << 53);
    }
    FUNC_ENABLE_IFSAME2(num_t,double,word_t,u64,num_t) inline randNum()
    {
        return (rng.next() >> 11) / (double)(1LL << 53);
    }
    inline void randGaussianPair(num_t& z1, num_t& z2)
    {
#if 0 // box muller transform
        num_t u1 = randNum();
        num_t u2 = (2.0*M_PI)*randNum();
        num_t r = sqrt(-2.0*log(u1));
        num_t s,c;
        sincosg(u2,&s,&c);
        z1 = r*c;
        z2 = r*s;
#else // marsaglia polar method
        num_t v1,v2,s,r;
        do
        {
            v1 = 2.0*randNum() - 1.0;
            v2 = 2.0*randNum() - 1.0;
            s = v1*v1 + v2*v2;
        }
        while (s >= 1.0);
        r = sqrt(-2.0*log(s)/s);
        z1 = v1*r;
        z2 = v2*r;
#endif
    }
    inline num_t randGaussian()
    {
        num_t z1,z2;
        randGaussianPair(z1,z2);
        return z1;
    }
    // random point in the biunit square/cube/hypercube
    inline point_t randPoint()
    {
        num_t x[dims];
        for (size_t i = 0; i < dims; ++i)
            x[i] = 2.0*randNum() - 1.0;
        return point_t(x);
    }
    // random point on the unit circle/sphere/hypersphere surface
    ENABLE_IFEQ(dims,1,point_t) inline randDirection()
    {
        static const num_t table[2] = {-1.0,1.0};
        return point_t(table[randBool()]);
    }
    ENABLE_IFEQ(dims,2,point_t) inline randDirection()
    {
        num_t a = (2.0*M_PI) * randNum();
        num_t sa,ca;
        sincosg(a,&sa,&ca);
        return point_t(ca,sa);
    }
    ENABLE_IFEQ(dims,3,point_t) inline randDirection()
    {
        num_t p = acos(2.0*randNum()-1.0);
        num_t t = (2.0*M_PI) * randNum();
        num_t x = sin(t)*cos(p);
        num_t y = sin(t)*sin(p);
        num_t z = cos(t);
        return point_t(x,y,z);
    }
    ENABLE_IF(dims>3,point_t) inline randDirection()
    {
        num_t x[dims];
        num_t dummy;
        for (size_t i = 0; i+1 < dims; i += 2)
            randGaussianPair(x[i],x[i+1]);
        if (dims % 2)
            x[dims-1] = randGaussian();
        point_t p(x);
        return p/p.norm2();
    }
    inline const Affine<num_t,2>& getPreAffine() const
    {
        return xf->getPreAffine();
    }
    inline const Affine<num_t,2>& getPostAffine() const
    {
        return xf->getPostAffine();
    }
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

template <typename num_t, size_t dims, typename rand_t>
struct XFormVar
{
    // function pointer
    std::function<void(IterState<num_t,dims,rand_t>&,const num_t*)> func;
    size_t index; // index of first variation parameter in varp (XForm class)
    // parameters are taken in order starting from index
    // the varp vector keeps the parameters compactly in memory
};

// xform (including final xform)
template <typename num_t, size_t dims, typename rand_t>
class XForm
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
                //vars2d<num_t,rand_t>::data.at(name);
                Variations<num_t,dims,rand_t>::get(name);
            var.func = varinfo.getFPtr();
            var.index = varp.size();
            num_t weight = varj["weight"].floatValue();
            if (weight == 0.0)
                continue; // skip zero weight variations
            vars.push_back(var);
            if (varinfo.getPPtr()) // store other parameters
                varinfo.getPPtr()(*this,varj,weight,varp);
            else // store just the weight by default
                varp.push_back(weight);
            pc_flags |= varinfo.getPCFlags();
        }
    }
    // optimize xform
    void optimize()
    {
        // zero weight variations are excluded in constructor
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
    inline const std::vector<XFormVar<num_t,dims,rand_t>>& getVariations() const
    {
        return vars;
    }
    inline const std::vector<num_t>& getVariationParams() const
    {
        return varp;
    }
    // iterate a state for the rendering process
    ENABLE_IF(dims<=2,void) inline applyIteration(
            IterState<num_t,dims,rand_t>& state) const
    {
        // for 2d, faster to apply identity affine than to branch
        state.xf = this;
        state.t = pre.apply_to(state.p);
        state.v = Point<num_t,dims>();
        for (XFormVar<num_t,dims,rand_t> v : vars)
            v.func(state,varp.data()+v.index);
        state.p = post.apply_to(state.v);
    }
    ENABLE_IF(dims>2,void) inline applyIteration(
            IterState<num_t,dims,rand_t>& state) const
    {
        // for 3d, probably faster to skip identity affine
        state.xf = this;
        if (has_pre)
            state.t = pre.apply_to(state.p);
        state.v = Point<num_t,dims>();
        for (XFormVar<num_t,dims,rand_t> v : vars)
            v.func(state,varp.data()+v.index);
        if (has_post)
            state.p = post.apply_to(state.v);
    }
};

// flame fractal
template <typename num_t, size_t dims, typename rand_t>
class Flame
{
private:
    std::array<size_t,dims> size;
    std::array<std::pair<num_t,num_t>,dims> bounds;
    std::vector<XForm<num_t,dims,rand_t>> xforms;
    XForm<num_t,dims,rand_t> final_xform;
    bool has_final_xform;
public:
    Flame(){}
    // construct from JSON data
    // throws error if something goes wrong
    Flame(const Json& input)
    {
        if (input["dimensions"].intValue() != dims)
            throw std::runtime_error("flame: dimension mismatch");
        JsonArray sizej = input["size"].arrayValue();
        JsonArray boundsj = input["bounds"].arrayValue();
        if (sizej.size() != dims)
            throw std::runtime_error("flame: incorrect size length");
        if (boundsj.size() != dims)
            throw std::runtime_error("flame: incorrect bounds length");
        num_t M = max_rect<num_t>::value;
        for (size_t i = 0; i < dims; ++i)
        {
            size[i] = sizej[i].floatValue();
            if (size[i] == 0 || size[i] > max_dim)
                throw std::runtime_error("flame: size not in [1,65535]");
            JsonArray boundj = boundsj[i].arrayValue();
            if (boundj.size() != 2)
                throw std::runtime_error("flame: incorrect bound format");
            num_t lo = boundj[0].floatValue();
            num_t hi = boundj[1].floatValue();
            bounds[i] = std::make_pair(lo,hi);
            if (lo < -M || lo > M || hi < -M || hi > M)
                throw std::runtime_error("flame: bound out of range");
            if (lo >= hi)
                throw std::runtime_error("flame: bound low >= bound high");
        }
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
    inline const std::array<size_t,dims>& getSize() const
    {
        return size;
    }
    inline const std::array<std::pair<num_t,num_t>,dims>& getBounds() const
    {
        return bounds;
    }
    inline const std::vector<XForm<num_t,dims,rand_t>>& getXForms() const
    {
        return xforms;
    }
    inline bool hasFinalXForm() const
    {
        return has_final_xform;
    }
    inline const XForm<num_t,dims,rand_t>& getFinalXForm() const
    {
        return final_xform;
    }
};

// render histogram only (count of samples in each pixel)
template <typename num_t, typename hist_t, typename cache_t, typename rand_t>
class RendererBasic
{
    static_assert(sizeof(hist_t) > sizeof(cache_t));
private:
    static constexpr hist_t cache_max = max_int_as<cache_t,hist_t>::value;
    Flame<num_t,2,rand_t> flame;
    hist_t *histogram;
    // main write location for renderer threads (smaller fits in cache better)
    // on overflow, add to histogram
    cache_t *bufcache;
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
    RendererBasic(const Flame<num_t,2,rand_t>& flame, hist_t *buf = nullptr):
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
        bufcache = new cache_t[X*Y]();
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
                __atomic_fetch_add(histogram+index,cache_max,
                    __ATOMIC_RELAXED);
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
    inline const Flame<num_t,2,rand_t>& getFlame() const
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
}

#undef FUNC_ENABLE_IF
#undef FUNC_ENABLE_IF2
#undef ENABLE_IF
#undef ENABLE_IFEQ
#undef likely
#undef unlikely
