#include "image.hpp"

namespace tkoz::flame
{

template <typename img_t>
bool writePng(const img_t& img, std::ostream& os)
{
    auto view = boost::gil::const_view(img);
    boost::gil::write_view(os,view,boost::gil::png_tag());
    return os.good();
}

// boost does not compile this with mono_img_t
//template bool writePng(const mono_img_t&,std::ostream&);
template bool writePng<gray_img_t<u8>>(const gray_img_t<u8>&,std::ostream&);
template bool writePng<gray_img_t<u16>>(const gray_img_t<u16>&,std::ostream&);
template bool writePng<rgb_img_t<u8>>(const rgb_img_t<u8>&,std::ostream&);
template bool writePng<rgb_img_t<u16>>(const rgb_img_t<u16>&,std::ostream&);

bool writePbm(const mono_img_t& img, std::ostream& os)
{
    size_t X = img.dimensions().x;
    size_t Y = img.dimensions().y;
    auto view = boost::gil::const_view(img);
    os << "P4\n" << X << " " << Y << "\n";
    for (size_t y = 0; y < Y; ++y)
    {
        auto row = view.row_begin(y);
        size_t x = 0;
        // write groups of 8 pixels as 1 byte each
        for (; x+7 < X; x += 8)
        {
            u8 c = 0;
            u8 b = (1 << 7);
            for (size_t i = x; i < x+8; ++i)
            {
                if (row[i] == 0)
                    c |= b;
                b >>= 1;
            }
            os.put(c);
        }
        // extra pixels
        if (x < X)
        {
            u8 c = 0;
            u8 b = (1 << 7);
            for (; x < X; ++x)
            {
                if (row[x] == 0)
                    c |= b;
                b >>= 1;
            }
            os.put(c);
        }
    }
    return os.good();
}

template <typename pix_t>
bool writePgm(const gray_img_t<pix_t>& img, std::ostream& os)
{
    static_assert(std::is_same_v<pix_t,u8> || std::is_same_v<pix_t,u16>);
    size_t X = img.dimensions().x;
    size_t Y = img.dimensions().y;
    auto view = boost::gil::const_view(img);
    os << "P5\n" << X << " " << Y << "\n";
    if (std::is_same_v<pix_t,u8>) // 8 bit
    {
        os << "255\n";
        for (size_t y = 0; y < Y; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x = 0; x < X; ++x)
                os.put(row[x]);
        }
    }
    else // 16 bit
    {
        os << "65535\n";
        for (size_t y = 0; y < Y; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x = 0; x < X; ++x)
            {
                os.put(row[x] >> 8);
                os.put(row[x]);
            }
        }
    }
    return os.good();
}

template bool writePgm<u8>(const gray_img_t<u8>&,std::ostream&);
template bool writePgm<u16>(const gray_img_t<u16>&,std::ostream&);

template <typename pix_t>
bool writePpm(const rgb_img_t<pix_t>& img, std::ostream& os)
{
    static_assert(std::is_same_v<pix_t,u8> || std::is_same_v<pix_t,u16>);
    size_t X = img.dimensions().x;
    size_t Y = img.dimensions().y;
    auto view = boost::gil::const_view(img);
    os << "P6\n" << X << " " << Y << "\n";
    if (std::is_same_v<pix_t,u8>) // 8 bit
    {
        os << "255\n";
        for (size_t y = 0; y < Y; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x = 0; x < X; ++x)
            {
                os.put(row[x][0]);
                os.put(row[x][1]);
                os.put(row[x][2]);
            }
        }
    }
    else // 16 bit
    {
        os << "65535\n";
        for (size_t y = 0; y < Y; ++y)
        {
            auto row = view.row_begin(y);
            for (size_t x  = 0; x < X; ++x)
            {
                os.put(row[x][0] >> 8);
                os.put(row[x][0]);
                os.put(row[x][1] >> 8);
                os.put(row[x][1]);
                os.put(row[x][2] >> 8);
                os.put(row[x][2]);
            }
        }
    }
    return os.good();
}

template bool writePpm<u8>(const rgb_img_t<u8>&,std::ostream&);
template bool writePpm<u16>(const rgb_img_t<u16>&,std::ostream&);

} // namespace tkoz::flame
