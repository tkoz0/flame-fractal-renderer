/*
SFINAE code
*/

#pragma once

#include <cstdlib>

// variadic template all same type
// https://www.fluentcpp.com/2019/01/25/variadic-number-function-parameters-type/
template <bool...> struct bool_pack{};
template <typename...Ts>
using vconj = std::is_same<bool_pack<true,Ts::value...>,
                           bool_pack<Ts::value...,true>>;
// all of Ts same as T
template <typename T, typename...Ts>
using vallsame = vconj<std::is_same<T,Ts>...>;
// all of Ts convertible to T
template <typename T, typename...Ts>
using vallconv = vconj<std::is_convertible<T,Ts>...>;

// similar to std::enable_if but for > 2 cases of function specialization
template <size_t N1, size_t N2, typename T = void>
struct enable_if_eq {};
template <size_t N, typename T>
struct enable_if_eq<N,N,T> { typedef T type; };
#define ENABLE_IFEQ(N1,N2,RET) template <typename RET2 = RET> \
    typename enable_if_eq<N1,N2,RET2>::type

