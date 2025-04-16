/*
Representation of affine transformation
*/

#pragma once

#include "point.hpp"

#include "../utils/json.hpp"

#include <array>
#include <cstdlib>
#include <iostream>

namespace tkoz::flame
{

// affine transformation in N dimensions
template <typename T, size_t N>
class Affine
{
    static_assert(N > 0 && N < 65536);

private:

    std::array<Point<T,N>,N> A; // linear transformation A*x + b
    Point<T,N> b;               // A as array of its rows

public:

    Affine(): A(), b() // default is identity transformation
    {
        // identity transformation
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
        if (!j.isObject())
            throw JsonError("Affine(Json&): not an object");
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
            if (!jo["A"].isArray())
                throw JsonError("Affine(Json&): A is not an array");
            JsonArray ja = jo["A"].arrayValue();
            if (ja.size() != N)
                throw JsonError("Affine(Json&): A is wrong size");
            for (size_t i = 0; i < N; ++i)
            {
                try
                {
                    A[i] = Point<T,N>(ja[i]);
                }
                catch (const std::exception& e)
                {
                    throw JsonError("Affine(Json&): error parsing A["
                        + std::to_string(i) + "]: " + e.what());
                }
            }
        }
        if (jo.find("b") == jo.end()) // use 0
            b = Point<T,N>();
        else
        {
            try
            {
                b = Point<T,N>(jo["b"]);
            }
            catch (const std::exception& e)
            {
                throw JsonError("Affine(Json&): error parsing b: "
                    + std::string(e.what()));
            }
        }
    }

    [[nodiscard]] inline const std::array<Point<T,N>,N>& getA() const
    {
        return A;
    }

    [[nodiscard]] inline const Point<T,N>& getB() const
    {
        return b;
    }

    [[nodiscard]] inline Point<T,N> apply_to(const Point<T,N>& x) const
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

} // namespace tkoz::flame
