#include "misc.hpp"

#include <cstdarg>

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
