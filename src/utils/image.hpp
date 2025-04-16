#pragma once

#include "../types/types.hpp"

#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>

#include <iostream>

namespace tkoz::flame
{

// boost monochrome image type
typedef boost::gil::gray1_image_t mono_img_t;

// boost gray image type
template <typename pix_t> struct gray_img;
template <> struct gray_img<u8> { typedef boost::gil::gray8_image_t type; };
template <> struct gray_img<u16> { typedef boost::gil::gray16_image_t type; };
template <typename pix_t> using gray_img_t = gray_img<pix_t>::type;

// boost rgb image type
template <typename pix_t> struct rgb_img;
template <> struct rgb_img<u8> { typedef boost::gil::rgb8_image_t type; };
template <> struct rgb_img<u16> { typedef boost::gil::rgb16_image_t type; };
template <typename pix_t> using rgb_img_t = rgb_img<pix_t>::type;

template <typename img_t>
bool writePng(const img_t& img, std::ostream& os);

bool writePbm(const mono_img_t& img, std::ostream& os);

template <typename pix_t>
bool writePgm(const gray_img_t<pix_t>& img, std::ostream& os);

template <typename pix_t>
bool writePpm(const rgb_img_t<pix_t>& img, std::ostream& os);

} // namespace tkoz::flame
