#include "utils.hpp"

#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>

// returns the nanosecond (or most precise) performance counter
size_t clock_nanotime()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return 1000000000uLL*t.tv_sec + t.tv_nsec;
}

// returns a string representing the entire contents of a file
std::string read_text_file(const std::string& name)
{
    std::ifstream input;
    input.open(name);
    std::stringstream ss;
    ss << input.rdbuf();
    input.close();
    return ss.str();
}

// print to stderr and exit with code 1
void error_and_exit(const char *f, ...)
{
    va_list args;
    va_start(args,f);
    vfprintf(stderr,f,args);
    va_end(args);
    exit(1);
}

// true if string ends with the provided suffix
bool string_ends_with(const std::string& string, const std::string& suffix)
{
    size_t l1 = string.length();
    size_t l2 = suffix.length();
    return l1 >= l2 && string.substr(l1-l2) == suffix;
}

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
