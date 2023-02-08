#include <iostream>

#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>

// write 8 bit pgm image
bool write_pgm(std::ostream& os, size_t X, size_t Y, uint8_t *img)
{
    os << "P5" << std::endl;
    os << X << " " << Y << std::endl;
    os << "255" << std::endl;
    os.write((char*)img,X*Y*sizeof(uint8_t));
    return os.good();
}

// write 16 bit pgm image
bool write_pgm(std::ostream& os, size_t X, size_t Y, uint16_t *img)
{
    os << "P5" << std::endl;
    os << X << " " << Y << std::endl;
    os << "65535" << std::endl;
    os.write((char*)img,X*Y*sizeof(uint16_t));
    return os.good();
}

// write 8 bit grayscale png image
bool write_png(std::ostream& os, size_t X, size_t Y, uint8_t *img)
{
    boost::gil::gray8_image_t image(X,Y);
    auto image_view = boost::gil::view(image);
    for (size_t y = 0; y < Y; ++y)
    {
        auto ptr = image_view.row_begin(y);
        for (size_t x = 0; x < X; ++x)
            ptr[x] = *(img++);
    }
    boost::gil::write_view(os,image_view,boost::gil::png_tag());
    return os.good();
}

// write 16 bit grayscale png image
bool write_png(std::ostream& os, size_t X, size_t Y, uint16_t *img)
{
    boost::gil::gray16_image_t image(X,Y);
    auto image_view = boost::gil::view(image);
    for (size_t y = 0; y < Y; ++y)
    {
        auto ptr = image_view.row_begin(y);
        for (size_t x = 0; x < X; ++x)
            ptr[x] = *(img++);
    }
    boost::gil::write_view(os,image_view,boost::gil::png_tag());
    return os.good();
}
