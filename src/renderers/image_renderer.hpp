#pragma once

#include "buffer_renderer.hpp"

#include "../utils/color.hpp"
#include "../utils/endian.hpp"
#include "../utils/image.hpp"
#include "../utils/scalers.hpp"

#include <boost/gil.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace tkoz::flame
{

/*
Image Renderer
- basic implementations so far
- supports 2d buffer only for now

TODO
- work out how to support 3d buffers (and 1d too)
- oversampling
- (variable radius) blur
- interpolation and super sampling
- brightness and vibrancy (see apophysis)
- per color channel adjustments
- alternative color spaces
- see apophysis for gamma threshold and filter radius
*/
template <bool enable_color = true>
class ImageRenderer
{
private:

    // point type
    typedef Point<num_t,2> point_t;
    // base variation type
    typedef vars::Variation<2> var_t;
    // xform type
    typedef XForm<2> xform_t;
    // flame type
    typedef Flame<2> flame_t;
    // number pair
    typedef std::pair<num_t,num_t> num_pair_t;

    // renderer object with buffer
    BufferRenderer<2,enable_color> buf_ren;

    // buffer viewer object for the filtered
    struct buf_view_t
    {
        const buf_elem_t *ptr;
        size_t xl,yl,r,cellsize;
        buf_view_t(const BufferRenderer<2,enable_color>& buf_ren)
        {
            ptr = buf_ren.getBuffer().data();
            xl = buf_ren.getFlame().getSize()[0];
            yl = buf_ren.getFlame().getSize()[1];
            r = buf_ren.getColorDims();
            cellsize = r+1;
        }
    };

public:

    ImageRenderer(const flame_t& flame): buf_ren(flame) {}

    ImageRenderer(const Json& json): buf_ren(json) {}

    // reference to buffer renderer used
    [[nodiscard]] inline BufferRenderer<2,enable_color>& getBufferRenderer()
    {
        return buf_ren;
    }

    // const reference to buffer renderer used
    [[nodiscard]] inline
    const BufferRenderer<2,enable_color>& getBufferRenderer() const
    {
        return buf_ren;
    }

    // buffer length in columns (X)
    [[nodiscard]] inline size_t getBufferSizeX() const
    {
        return buf_ren.getFlame().getSize()[0];
    }

    // buffer length in rows (Y)
    [[nodiscard]] inline size_t getBufferSizeY() const
    {
        return buf_ren.getFlame().getSize()[1];
    }

    // buffer size (both dimensions)
    [[nodiscard]] inline std::pair<size_t,size_t> getBufferSize() const
    {
        return std::make_pair(getBufferSizeX(),getBufferSizeY());
    }

    // number of color dimensions
    [[nodiscard]] inline size_t getColorDims() const
    {
        return buf_ren.getColorDims();
    }

    // get min and max value of a function over all pixels
    template <typename T>
    [[nodiscard]] std::pair<T,T> getValueBounds(
        util::cell_func_t<T> func) const
    {
        buf_view_t v(buf_ren);
        T min,max;
        min = max = func(v.ptr,v.r);
        for (size_t i = 0; i < buf_ren.getBufferNumCells(); ++i)
        {
            T val = func(v.ptr,v.r);
            min = std::min(min,val);
            max = std::max(max,val);
            v.ptr += v.cellsize;
        }
        return std::make_pair(min,max);
    }

    // render binary image (monochrome bitmap)
    // bool function determines which points are white
    [[nodiscard]] mono_img_t renderBinaryImage(
        util::cell_func_t<bool> white) const
    {
        buf_view_t v(buf_ren);
        mono_img_t ret(v.xl,v.yl);
        auto view = boost::gil::view(ret);
        for (size_t y = 0; y < v.yl; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x = 0; x < v.xl; ++x)
            {
                row[x] = white(v.ptr,v.r) ? 255 : 0;
                v.ptr += v.cellsize;
            }
        }
        return ret;
    }

    // render gray image
    template <typename pix_t> [[nodiscard]]
    gray_img_t<pix_t> renderGrayImage(
        util::cell_func_t<num_t> scaler) const
    {
        buf_view_t v(buf_ren);
        gray_img_t<pix_t> ret(v.xl,v.yl);
        auto view = boost::gil::view(ret);
        for (size_t y = 0; y < v.yl; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x = 0; x < v.xl; ++x)
            {
                num_t l = scaler(v.ptr,v.r);
                row[x] = (pix_t)(l * pix_scale_v<pix_t,num_t>);
                v.ptr += v.cellsize;
            }
        }
        return ret;
    }

    // render color image
    template <typename pix_t> [[nodiscard]]
    rgb_img_t<pix_t> renderColorImageRGB(
            util::cell_func_t<std::tuple<num_t,num_t,num_t>> colorer) const
    {
        buf_view_t v(buf_ren);
        rgb_img_t<pix_t> ret(v.xl,v.yl);
        auto view = boost::gil::view(ret);
        for (size_t y = 0; y < v.yl; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x = 0; x < v.xl; ++x)
            {
                auto [r,g,b] = colorer(v.ptr,v.r);
                pix_t rp = (pix_t)(r * pix_scale_v<pix_t,num_t>);
                pix_t gp = (pix_t)(g * pix_scale_v<pix_t,num_t>);
                pix_t bp = (pix_t)(b * pix_scale_v<pix_t,num_t>);
                row[x] = {rp,gp,bp};
                v.ptr += v.cellsize;
            }
        }
        return ret;
    }
};

} // namespace tkoz::flame
