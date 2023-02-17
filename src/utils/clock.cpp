#include <cstdlib>
#include <ctime>

#include "clock.hpp"

// returns the nanosecond (or most precise) performance counter
size_t clock_nanotime()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return 1000000000uLL*t.tv_sec + t.tv_nsec;
}
