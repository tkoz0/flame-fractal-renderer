/*
SFINAE code
for some complicated stuff with advanced features
*/

#pragma once

#include <cstdlib>
#include <type_traits>

// variadic template all same type
// https://www.fluentcpp.com/2019/01/25/variadic-number-function-parameters-type/
template <bool...> struct bool_pack{};
template <typename...Ts>
using vconj = std::is_same<bool_pack<true,Ts::value...>,
                           bool_pack<Ts::value...,true>>;
template <typename...Ts>
static constexpr auto vconj_v = vconj<Ts...>::value;

// all of Ts same as T
template <typename T, typename...Ts>
using vallsame = vconj<std::is_same<T,Ts>...>;
template <typename T, typename...Ts>
static constexpr auto vallsame_v = vallsame<T,Ts...>::value;

// all of Ts convertible to T
template <typename T, typename...Ts>
using vallconv = vconj<std::is_convertible<T,Ts>...>;
template <typename T, typename...Ts>
static constexpr auto vallconv_v = vallconv<T,Ts...>::value;

// similar to std::enable_if but for > 2 cases of function specialization
template <size_t N1, size_t N2, typename T = void>
struct enable_if_eq {};
template <size_t N, typename T>
struct enable_if_eq<N,N,T> { typedef T type; };
template <size_t N1, size_t N2, typename T = void>
using enable_if_eq_t = enable_if_eq<N1,N2,T>::type;
