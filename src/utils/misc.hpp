#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// returns a string representing the entire contents of a file
std::string read_text_file(const std::string& name);

// print to stderr and exit with code 1
void error_and_exit(const char *f, ...);

// true if string ends with the provided suffix
bool string_ends_with(const std::string& string, const std::string& suffix);
