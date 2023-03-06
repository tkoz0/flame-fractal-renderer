#pragma once

#include <iostream>

#include "endian.hpp"

// write 8 bit pgm image
bool write_pgm(std::ostream& os, size_t X, size_t Y, uint8_t *img);

// write 16 bit pgm image
bool write_pgm(std::ostream& os, size_t X, size_t Y, uint16_t *img);

// write 8 bit grayscale png image
bool write_png(std::ostream& os, size_t X, size_t Y, uint8_t *img);

// write 16 bit grayscale png image
bool write_png(std::ostream& os, size_t X, size_t Y, uint16_t *img);
