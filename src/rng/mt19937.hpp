/*
MT19937 Mersenne twister implementation
*/

#pragma once

#include <array>
#include <cstdint>
#include <ctime>

#include "../utils/clock.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define FUNC_ENABLE_IF(T1,T2,RET) template <typename dummy = T1> \
    typename std::enable_if<std::is_same<dummy,T2>::value,RET>::type

// integer types
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;

template <typename T> struct mt_params {};
template <> struct mt_params<u32>
{
    static const size_t w = 32;
    static const size_t n = 624;
    static const size_t m = 397;
    static const size_t r = 31;
    static const u32 a = 0x9908B0DF;
    static const u32 u = 11;
    static const u32 d = 0xFFFFFFFF;
    static const u32 s = 7;
    static const u32 b = 0x9D2C5680;
    static const u32 t = 15;
    static const u32 c = 0xEFC60000;
    static const u32 l = 18;
    static const u32 f = 1812433253;
    static const u32 seed_uniq = 1853726117;
};

template <> struct mt_params<u64>
{
    static const size_t w = 64;
    static const size_t n = 312;
    static const size_t m = 156;
    static const size_t r = 31;
    static const u64 a = 0xB5026F5AA96619E9;
    static const u64 u = 29;
    static const u64 d = 0x5555555555555555;
    static const u64 s = 17;
    static const u64 b = 0x71D67FFFEDA60000;
    static const u64 t = 37;
    static const u64 c = 0xFFF7EEE000000000;
    static const u64 l = 43;
    static const u64 f = 6364136223846793005;
    static const u64 seed_uniq = 7460189664030870853;
};

template <typename T> class MT19937
{
private:
    typedef mt_params<T> params;
    std::array<T,params::n> MT;
    size_t index;
    static const T mask_lo = (((T)1) << params::r) - 1;
    static const T mask_hi = ~mask_lo;
    // seed must be nonzero
    inline void _seed(T seed)
    {
        index = params::n;
        MT[0] = seed;
        for (size_t i = 1; i < params::n; ++i)
            MT[i] = params::f * (MT[i-1] ^ (MT[i-1] >> (params::w - 2))) + i;
    }
    inline void _twist()
    {
        for (size_t i = 0; i < params::n; ++i)
        {
            T x = (MT[i] & mask_hi) | (MT[(i+1) % params::n] & mask_lo);
            T xA = (x >> 1) ^ (params::a * (x & 1));
            MT[i] = MT[(i+params::m) % params::n] ^ xA;
        }
        index = 0;
    }
    inline T _get()
    {
        if (unlikely(index >= params::n))
            _twist();
        T y = MT[index++];
        y ^= (y >> params::u) & params::d;
        y ^= (y << params::s) & params::b;
        y ^= (y << params::t) & params::c;
        y ^= y >> params::l;
        return y;
    }
public:
    inline MT19937(T seed = 0)
    {
        setSeed(seed);
    }
    inline void setSeed(T seed = 0)
    {
        static T seed_uniq = params::seed_uniq;
        if (seed)
            _seed(seed);
        else
        {
            seed_uniq = (seed_uniq * 711671537) + 1002671653;
            _seed(clock_nanotime() ^ seed_uniq ^ time(NULL));
        }
    }
    inline T nextWord()
    {
        return _get();
    }
    FUNC_ENABLE_IF(T,u32,float) inline nextFloat()
    {
        return (_get() >> 8) / (float) 0x1000000;
    }
    FUNC_ENABLE_IF(T,u64,double) inline nextFloat()
    {
        return (_get() >> 11) / (double) 0x20000000000000L;
    }
};

extern template class MT19937<u32>;
extern template class MT19937<u64>;

#undef likely
#undef unlikely
#undef FUNC_ENABLE_IF
