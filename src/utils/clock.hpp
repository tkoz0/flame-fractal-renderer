#pragma once

#include <cstdlib>
#include <ctime>

namespace tkoz::flame::util
{

// returns the nanosecond (or most precise) performance counter
[[nodiscard]] inline size_t clock_nanotime()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return 1000000000uLL*t.tv_sec + t.tv_nsec;
}

} // namespace tkoz::flame::util
