#pragma once

#include <cstdlib>
#include <thread>

// get number of threads supported
inline size_t number_of_threads()
{
    size_t n = std::thread::hardware_concurrency();
    return std::max((size_t)1,n); // return 1 if n is 0
}
