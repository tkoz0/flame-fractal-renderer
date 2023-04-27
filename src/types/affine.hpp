/*
Representation of affine transformation
*/

#pragma once

#include <array>
#include <cstdlib>
#include <iostream>

// forward declaration
namespace tkoz::flame
{
template <typename T, size_t N> class Affine;
}

#include "../utils/json_small.hpp"
#include "point.hpp"

namespace tkoz::flame
{

// affine transformation in N dimensions
template <typename T, size_t N>
class Affine
{
    static_assert(N > 0 && N < 100);
private:
    std::array<Point<T,N>,N> A; // linear transformation A*x + b
    Point<T,N> b;               // A as array of its rows
public:
    Affine(): A(), b()
    {
        for (size_t i = 0; i < N; ++i)
            A[i][i] = 1;
    }
    Affine(const T *A_, const T *b_) // from length N*N and N arrays
    {
        b = Point<T,N>(b_);
        for (size_t i = 0; i < N; ++i)
            A[i] = Point<T,N>(A_+(N*i));
    }
    Affine(const Json& j) // {"A": [2d array], "b": [1d array]}}
    {
        JsonObject jo = j.objectValue();
        if (jo.find("A") == jo.end()) // use identity matrix
        {
            for (size_t i = 0; i < N; ++i)
            {
                A[i] = Point<T,N>();
                A[i][i] = 1;
            }
        }
        else
        {
            JsonArray ja = jo["A"].arrayValue();
            if (ja.size() != N)
                throw std::runtime_error("affine: incorrect size");
            for (size_t i = 0; i < N; ++i)
                A[i] = Point<T,N>(ja[i]);
        }
        if (jo.find("b") == jo.end()) // use 0
            b = Point<T,N>();
        else
            b = Point<T,N>(jo["b"]);
    }
    inline const std::array<Point<T,N>,N>& getA() const
    {
        return A;
    }
    inline const Point<T,N>& getB() const
    {
        return b;
    }
    inline Point<T,N> apply_to(const Point<T,N>& x) const
    {
        Point<T,N> ret(b);
        for (size_t i = 0; i < N; ++i)
            ret[i] += A[i] * x;
        return ret;
    }
};

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Affine<T,N>& a)
{
    os << "{A=[" << a[0];
    for (size_t i = 1; i < N; ++i)
        os << "," << a[i];
    os << "],b=" << a.getB() << "}";
    return os;
}

}
