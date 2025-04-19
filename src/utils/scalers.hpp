#pragma once

#include "../types/types.hpp"

#include <functional>

namespace tkoz::flame::util
{

// function type for operating on a buffer cell
template <typename ret_t>
using cell_func_t = std::function<ret_t(const buf_elem_t*,size_t r)>;

// scaler function type, map hist_t to num_t
// must be nondecreasing and have range [0, max num_t value]
typedef std::function<num_t(hist_t)> scaler_t;

scaler_t scale_log = [](hist_t h) -> num_t
{
    return std::log(1.0 + (num_t)h);
};

scaler_t scale_linear = [](hist_t h) -> num_t
{
    return (num_t)h;
};

scaler_t scale_binary = [](hist_t h) -> num_t
{
    return h > 0 ? 1.0 : 0.0;
};

template <u32 root>
scaler_t scale_root = [](hist_t h) -> num_t
{
    static_assert(root > 0);
    if (root == 1)
        return (num_t)h;
    if (root == 2)
        return std::sqrt(h);
    if (root == 3)
        return std::cbrt(h);
    return std::pow(h,1.0/root);
};

/*
further function ideas (p positive integer, c positive real)

x -> log(1 + x)
x -> x
x -> 0 if x = 0, 1 if x > 0
x -> 0 if x <= c, 1 if x > c
x -> x^(1/p)
x -> x^p / (x^p + c)
x -> arctan(x^(1/p)/c)
x -> 1 - 1/(x + c)
x -> log(1 + x/c)
x -> log(1 + log(1 + x))
x -> log(1 + log(1 + x/c1)/c2)
x -> log(1 + x)^p
x -> log(1 + x/c)^p
x -> 1 - c^x where c in [0,1)
x -> x + c
x -> x^(1/p) + c

maybe allow functions to be decreasing or not even monotonic
*/

} // namespace tkoz::flame::util
