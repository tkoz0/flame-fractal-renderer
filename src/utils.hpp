#pragma once

#include <string>

// returns the nanosecond (or most precise) performance counter
size_t clock_nanotime();

// returns a string representing the entire contents of a file
std::string read_text_file(const std::string& name);

// print to stderr and exit with code 1
void error_and_exit(const char *f, ...);

// true if string ends with the provided suffix
bool string_ends_with(const std::string& string, const std::string& suffix);

// write 8 bit pgm image
bool write_pgm(std::ostream& os, size_t X, size_t Y, uint8_t *img);

// write 16 bit pgm image
bool write_pgm(std::ostream& os, size_t X, size_t Y, uint16_t *img);

// write 8 bit grayscale png image
bool write_png(std::ostream& os, size_t X, size_t Y, uint8_t *img);

// write 16 bit grayscale png image
bool write_png(std::ostream& os, size_t X, size_t Y, uint16_t *img);
