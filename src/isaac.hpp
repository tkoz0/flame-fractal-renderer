#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <tuple>
#include <type_traits>

#include "utils.hpp"

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

// macros for using SFINAE
#define FUNC_ENABLE_IF(T1,T2,RET) template <typename dummy = T1> \
    typename std::enable_if<std::is_same<dummy,T2>::value,RET>::type
//#define VAR_ENABLE_IF(T1,T2,TYPE) template <typename dummy = T1>
//static const typename std::enable_if<std::is_same<dummy,T2>::value,TYPE>::type

// integer types
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;

// golden ratio constant based on word type
template <typename T> struct golden_ratio {};
template <> struct golden_ratio<u8>
{ static const u8 value  = 0x9e; };
template <> struct golden_ratio<u16>
{ static const u16 value = 0x9e37; };
template <> struct golden_ratio<u32>
{ static const u32 value = 0x9e3779b9u; };
template <> struct golden_ratio<u64>
{ static const u64 value = 0x9e3779b97f4a7c13uLL; };

/*
ISAAC random number generator. Based on the code from
https://www.burtleburtle.net/bob/rand/isaacafa.html
word_t - word size (u32 or u64)
rparam - determines internal state size (2^rparam words)
*/
template <typename word_t, size_t rparam>
class Isaac
{
    static_assert(rparam >= 3 && rparam <= 16, "rparam out of range");
    static_assert(std::is_same<word_t,u32>::value
            || std::is_same<word_t,u64>::value, "word_t must be u32 or u64");
public:
    // number of words in the internal state
    static const size_t randsiz = 1 << rparam;
    static const size_t randsizmask = randsiz-1;
private:
    size_t randcnt; // number of unused values in randrsl
    word_t randrsl[randsiz]; // generated values
    word_t randmem[randsiz]; // internal state
    word_t randa,randb,randc;
    // advances internal state (generates 2^rparam new words)
    void gen()
    {
        word_t a,b,x,y,*m,*mm,*m2,*r,*mend;
        mm = randmem;
        r = randrsl;
        a = randa;
        b = randb + (++randc);
        for (m = mm, mend = m2 = m+(randsiz/2); m < mend;)
            rngstep4(a,b,mm,m,m2,r,x,y);
        for (m2 = mm; m2 < mend;)
            rngstep4(a,b,mm,m,m2,r,x,y);
        randa = a;
        randb = b;
    }
    // initialize state, flag = use randrsl,randa,randb,randc as seed
    void init(bool flag, word_t a0, word_t b0, word_t c0)
    {
        size_t i;
        word_t a,b,c,d,e,f,g,h;
        word_t *m = randmem;
        word_t *r = randrsl;
        randa = a0;
        randb = b0;
        randc = c0;
        a = b = c = d = e = f = g = h = golden_ratio<word_t>::value;
        mix(a,b,c,d,e,f,g,h); // scramble (4 times)
        mix(a,b,c,d,e,f,g,h);
        mix(a,b,c,d,e,f,g,h);
        mix(a,b,c,d,e,f,g,h);
        for (i = 0; i < randsiz; i += 8)
        {
            if (flag)
            {
                a += r[i+0]; b += r[i+1]; c += r[i+2]; d += r[i+3];
                e += r[i+4]; f += r[i+5]; g += r[i+6]; h += r[i+7];
            }
            mix(a,b,c,d,e,f,g,h);
            m[i+0] = a; m[i+1] = b; m[i+2] = c; m[i+3] = d;
            m[i+4] = e; m[i+5] = f; m[i+6] = g; m[i+7] = h;
        }
        if (flag) // 2nd pass to make seed affect all of randmem
        {
            for (i = 0; i < randsiz; i += 8)
            {
                a += m[i+0]; b += m[i+1]; c += m[i+2]; d += m[i+3];
                e += m[i+4]; f += m[i+5]; g += m[i+6]; h += m[i+7];
                mix(a,b,c,d,e,f,g,h);
                m[i+0] = a; m[i+1] = b; m[i+2] = c; m[i+3] = d;
                m[i+4] = e; m[i+5] = f; m[i+6] = g; m[i+7] = h;
            }
        }
        gen(); // initial results
        randcnt = randsiz;
    }
    // ind macro from the C code
    //static inline word_t ind(word_t *mm, word_t x);
    FUNC_ENABLE_IF(word_t,u32,u32) static inline ind(u32 *mm, u32 x)
    {
        return mm[(x >> 2) & randsizmask];
    }
    FUNC_ENABLE_IF(word_t,u64,u64) static inline ind(u64 *mm, u64 x)
    {
        return mm[(x >> 3) & randsizmask];
    }
    // rngstep macro from the C code
    static inline void rngstep(word_t mix, word_t& a, word_t& b, word_t *mm,
        word_t*& m, word_t*& m2, word_t*& r, word_t& x, word_t& y)
    {
        x = *m;
        a = mix + *(m2++);
        *(m++) = y = ind(mm,x) + a + b;
        *(r++) = b = ind(mm, y >> rparam) + x;
    }
    // mix macro from the C code
    //static inline void mix(word_t& a, word_t& b, word_t& c,
    //    word_t& d, word_t& e, word_t& f, word_t& g, word_t& h);
    FUNC_ENABLE_IF(word_t,u32,void) static inline mix(u32& a, u32& b,
        u32& c, u32& d, u32& e, u32& f, u32& g, u32& h)
    {
        a ^= b << 11; d += a; b += c;
        b ^= c >>  2; e += b; c += d;
        c ^= d <<  8; f += c; d += e;
        d ^= e >> 16; g += d; e += f;
        e ^= f << 10; h += e; f += g;
        f ^= g >>  4; a += f; g += h;
        g ^= h <<  8; b += g; h += a;
        h ^= a >>  9; c += h; a += b;
    }
    FUNC_ENABLE_IF(word_t,u64,void) static inline mix(u64& a, u64& b,
        u64& c, u64& d, u64& e, u64& f, u64& g, u64& h)
    {
        a -= e; f ^= h >>  9; h += a;
        b -= f; g ^= a <<  9; a += b;
        c -= g; h ^= b >> 23; b += c;
        d -= h; a ^= c << 15; c += d;
        e -= a; b ^= d >> 14; d += e;
        f -= b; c ^= e << 20; e += f;
        g -= c; d ^= f >> 17; f += g;
        h -= d; e ^= g << 14; g += h;
    }
    // combine the 4 rngstep uses for the generator loop
    //static inline void rngstep4(word_t& a, word_t& b, word_t *mm,
    //    word_t*& m, word_t*& m2, word_t*& r, word_t& x, word_t& y);
    FUNC_ENABLE_IF(word_t,u32,void) static inline rngstep4(u32& a, u32& b,
        u32 *mm, u32*& m, u32*& m2, u32*& r, u32& x, u32& y)
    {
        rngstep(a^(a<<13),a,b,mm,m,m2,r,x,y);
        rngstep(a^(a>> 6),a,b,mm,m,m2,r,x,y);
        rngstep(a^(a<< 2),a,b,mm,m,m2,r,x,y);
        rngstep(a^(a>>16),a,b,mm,m,m2,r,x,y);
    }
    FUNC_ENABLE_IF(word_t,u64,void) static inline rngstep4(u64& a, u64& b,
        u64 *mm, u64*& m, u64*& m2, u64*& r, u64& x, u64& y)
    {
        rngstep(~(a^(a<<21)),a,b,mm,m,m2,r,x,y);
        rngstep(  a^(a>> 5) ,a,b,mm,m,m2,r,x,y);
        rngstep(  a^(a<<12) ,a,b,mm,m,m2,r,x,y);
        rngstep(  a^(a>>33) ,a,b,mm,m,m2,r,x,y);
    }
    // copy constructor helper
    void copy(const Isaac& a)
    {
        randcnt = a.randcnt;
        memcpy(randrsl,a.randrsl,sizeof(word_t)*randsiz);
        memcpy(randmem,a.randmem,sizeof(word_t)*randsiz);
        randa = a.randa;
        randb = a.randb;
        randc = a.randc;
    }
public:
    // random seed based on system time
    Isaac()
    {
        setSeed();
    }
    // seed with a single 32 bit integer
    Isaac(u32 seed)
    {
        setSeed(seed);
    }
    // seed with a single 64 bit integer
    Isaac(u64 seed)
    {
        setSeed(seed);
    }
    // functions for single integer seed, dependent on word_t
    FUNC_ENABLE_IF(word_t,u32,void) setSeed()
    {
        static u32 seed_uniq = 2451404387u;
        seed_uniq = (seed_uniq * 229) + 137;
        setSeed(clock_nanotime() ^ seed_uniq ^ time(NULL));
    }
    FUNC_ENABLE_IF(word_t,u64,void) setSeed()
    {
        static u64 seed_uniq = 11713835213681433683uLL;
        seed_uniq = (seed_uniq * 53161) + 46457;
        setSeed(clock_nanotime() ^ seed_uniq ^ time(NULL));
    }
    FUNC_ENABLE_IF(word_t,u32,void) setSeed(u32 s)
    {
        setSeed(s,~s,(s<<16)|(s>>16));
    }
    FUNC_ENABLE_IF(word_t,u32,void) setSeed(u64 s)
    {
        setSeed(s,s>>32,s^(s>>32));
    }
    FUNC_ENABLE_IF(word_t,u64,void) setSeed(u32 s)
    {
        // 32 bit prime with 16 one bits
        setSeed((((u64)s)<<32)|(s^2451404387u));
    }
    FUNC_ENABLE_IF(word_t,u64,void) setSeed(u64 s)
    {
        // 64 bit prime with 32 one bits
        setSeed(s,~s,s^11713835213681433683uLL);
    }
    // main seed setting function
    void setSeed(word_t a0, word_t b0, word_t c0, word_t *seed = nullptr,
            bool force_flag = false)
    {
        if (seed)
            memcpy(randrsl,seed,sizeof(word_t)*randsiz);
        else
            memset(randrsl,0,sizeof(word_t)*randsiz);
        init(seed != nullptr || force_flag, a0, b0, c0);
    }
    // copy constructor
    Isaac(const Isaac& a)
    {
        copy(a);
    }
    // assignment operator
    Isaac& operator=(const Isaac& a)
    {
        if (&a != this)
            copy(a);
        return *this;
    }
    /*
    initializes the RNG state using 3 initial words and optionally an array
    if specified, `seed` array must point to `randsiz` words of data
    force_flag = use array of zeroes when `seed == nullptr`
    */
    Isaac(word_t a0, word_t b0, word_t c0,
        word_t *seed = nullptr, bool force_flag = false)
    {
        setSeed(a0,b0,c0,seed,force_flag);
    }
    /*
    initialize RNG state using only an array
    if `seed == nullptr`, an array is not used as a seed
    if specified, `seed` array must point to `randsiz` words of data
    this works the same as the other constructor but a0,b0,c0 are 0
    */
    Isaac(word_t *seed, bool force_flag = false)
    {
        setSeed(0,0,0,seed,force_flag);
    }
    // extract the next word, advancing the internal state if necessary
    inline word_t next()
    {
        if (unlikely(randcnt-- == 0))
        {
            gen();
            randcnt = randsiz-1;
        }
        return randrsl[randcnt];
    }
    // number of words left in the current internal state
    inline size_t count() const
    {
        return randcnt;
    }
    // advance the internal state (generates `randsiz` new words)
    inline void generate()
    {
        gen();
        randcnt = randsiz;
    }
    // returns the values array (length = randsiz)
    inline const word_t *valuesArray() const
    {
        return randrsl;
    }
    // returns the internal state array (length = randsiz)
    // random numbers should not be taken from this
    inline const word_t *stateArray() const
    {
        return randmem;
    }
    // returns the 3 internal state words
    inline std::tuple<word_t,word_t,word_t> stateWords() const
    {
        return std::make_tuple(randa,randb,randc);
    }
};

#undef FUNC_ENABLE_IF
#undef likely
#undef unlikely
