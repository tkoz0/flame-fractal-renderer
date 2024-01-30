#pragma once

#include <cmath>
#include <ctgmath>

#include "../types/constants.hpp"
#include "../utils/math.hpp"

namespace tkoz::flame
{

// forward declarations
struct color_lrgb;
struct color_srgb;
struct color_hsi;
struct color_hsl;
struct color_hsv;
struct color_xyz;
struct color_oklab;

// linear RGB
// r,g,b in [0,1]
struct color_lrgb
{
    num_t r,g,b;
    inline color_lrgb(num_t r, num_t g, num_t b): r(r), g(g), b(b) {}
    inline color_lrgb(color_hsl);
    inline color_lrgb(color_hsv);
    inline color_lrgb(color_hsi);
    inline color_lrgb(color_xyz);
    inline color_lrgb(color_srgb);
    inline color_lrgb(color_oklab);
};

// standard RGB (linear with gamma correction)
// r,g,b in [0,1]
struct color_srgb
{
    num_t r,g,b;
    inline color_srgb(num_t r, num_t g, num_t b): r(r), g(g), b(b) {}
    inline color_srgb(color_lrgb);
};

// h in [0,2pi) s,i in [0,1]
struct color_hsi
{
    num_t h,s,i;
    inline color_hsi(num_t h, num_t s, num_t i): h(h), s(s), i(i) {}
    inline color_hsi(color_lrgb);
};

// h in [0,2pi) s,l in [0,1]
struct color_hsl
{
    num_t h,s,l;
    inline color_hsl(num_t h, num_t s, num_t l): h(h), s(s), l(l) {}
    inline color_hsl(color_hsv);
    inline color_hsl(color_lrgb);
};

// h in [0,2pi) s,v in [0,1]
struct color_hsv
{
    num_t h,s,v;
    inline color_hsv(num_t h, num_t s, num_t v): h(h), s(s), v(v) {}
    inline color_hsv(color_hsl);
    inline color_hsv(color_lrgb);
};

// cie xyz x,y,z in [0,1]
struct color_xyz
{
    num_t x,y,z;
    inline color_xyz(num_t x, num_t y, num_t z): x(x), y(y), z(z) {}
    inline color_xyz(color_lrgb);
    inline color_xyz(color_oklab);
};

// oklab L,a,b in [0,1]
// https://bottosson.github.io/posts/oklab/
struct color_oklab
{
    num_t L,a,b;
    inline color_oklab(num_t L, num_t a, num_t b): L(L), a(a), b(b) {}
    inline color_oklab(color_xyz);
    inline color_oklab(color_lrgb);
};

const num_t hp_scaler = c_pi_3;

// helper for intermediate calculation to convert hsl/hsv to rgb
inline void _set_from_hue(num_t hp, num_t x, num_t c, num_t& r1, num_t& g1, num_t& b1)
{
    switch((u32)hp)
    {
    case 0: r1 = c; g1 = x; b1 = 0; break;
    case 1: r1 = x; g1 = c; b1 = 0; break;
    case 2: r1 = 0; g1 = c; b1 = x; break;
    case 3: r1 = 0; g1 = x; b1 = c; break;
    case 4: r1 = x; g1 = 0; b1 = c; break;
    case 5: r1 = c; g1 = 0; b1 = x; break;
    default:
        throw std::runtime_error("_set_from_hue: " + std::to_string(hp));
    }
}

// helper for intermediate calculation to convert rgb to hue
inline num_t _hue_from_rgb(num_t c, num_t v, num_t r, num_t g, num_t b)
{
    if (c == 0.0)
        return 0.0;
    if (v == r)
        return c_pi_3 * fmod((g-b)/c,6.0);
    if (v == g)
        return c_pi_3 * ((b-r)/c + 2.0);
    if (v == b)
        return c_pi_3 * ((r-g)/c + 4.0);
    throw std::runtime_error("_hue_from_rgb");
}

inline num_t _srgb_to_linear(num_t c)
{
    return c <= 0.04045 ? c/12.92 : std::pow((c+0.055)/1.055,2.4);
}

inline num_t _linear_to_srgb(num_t c)
{
    return c <= 0.0031308 ? 12.92*c : 1.055*std::pow(c,1.0/2.4) - 0.055;
}

inline color_lrgb::color_lrgb(color_hsl z)
{
    num_t c = (1.0 - std::abs(2*z.l - 1.0)) * z.s;
    num_t hp = z.h / hp_scaler;
    num_t x = c * (1.0 - std::abs(std::fmod(hp,2.0) - 1.0));
    num_t r1,g1,b1;
    _set_from_hue(hp,x,c,r1,g1,b1);
    num_t m = z.l - c/2.0;
    r = r1 + m;
    g = g1 + m;
    b = b1 + m;
}

inline color_lrgb::color_lrgb(color_hsv z)
{
    num_t c = z.v * z.s;
    num_t hp = z.h / hp_scaler;
    num_t x = c * (1.0 - std::abs(std::fmod(hp,2.0) - 1.0));
    num_t r1,g1,b1;
    _set_from_hue(hp,x,c,r1,g1,b1);
    num_t m = z.v - c;
    r = r1 + m;
    g = g1 + m;
    b = b1 + m;
}

inline color_lrgb::color_lrgb(color_hsi y)
{
    num_t hp = y.h / hp_scaler;
    num_t z = 1.0 - std::abs(std::fmod(hp,2.0) - 1.0);
    num_t c = (3.0 * y.i * y.s) / (1.0 + z);
    num_t x = c * z;
    num_t r1,g1,b1;
    _set_from_hue(hp,x,c,r1,g1,b1);
    num_t m = y.i * (1.0 - y.s);
    r = r1 + m;
    g = g1 + m;
    b = b1 + m;
}

inline color_lrgb::color_lrgb(color_xyz w)
{
    r = +3.2406255*w.x - 1.5372080*w.y - 0.4986286*w.z;
    g = -0.9689307*w.x + 1.8757561*w.y + 0.0415175*w.z;
    b = +0.0577101*w.x - 0.2040211*w.y + 1.0569959*w.z;
}

inline color_lrgb::color_lrgb(color_srgb z)
{
    r = _linear_to_srgb(z.r);
    g = _linear_to_srgb(z.g);
    b = _linear_to_srgb(z.b);
}

inline color_lrgb::color_lrgb(color_oklab z)
{
    num_t l = z.L + 0.3963377774*z.a + 0.2158037573*z.b;
    num_t m = z.L - 0.1055613458*z.a - 0.0638541728*z.b;
    num_t s = z.L - 0.0894841775*z.a - 1.2914855480*z.b;
    l = l*l*l;
    m = m*m*m;
    s = s*s*s;
    r = +4.0767416621*l - 3.3077115913*m + 0.2309699292*s;
    g = -1.2684380046*l + 2.6097574011*m - 0.3413193965*s;
    b = -0.0041960863*l - 0.7034186147*m + 1.7076147010*s;
}

inline color_srgb::color_srgb(color_lrgb z)
{
    r = _srgb_to_linear(z.r);
    g = _srgb_to_linear(z.g);
    b = _srgb_to_linear(z.b);
}

inline color_hsi::color_hsi(color_lrgb z)
{
    num_t rgbsum = z.r + z.g + z.b;
    i = rgbsum / 3.0;
    s = 1.0 - (3.0/rgbsum)*std::min(std::min(z.r,z.g),z.b);
    num_t rmg = z.r - z.g;
    num_t rmb = z.r - z.b;
    h = std::acos((rmg +(rmb))
        / (2.0 * std::sqrt(rmg*rmg + rmb*(z.g-z.b))));
    if (z.b > z.g)
        h = c_2pi - h;
}

inline color_hsl::color_hsl(color_hsv z)
{
    h = z.h;
    l = z.v * (1.0 - z.s/2.0);
    s = (l == 0.0 || l == 1.0) ? 0.0 : (z.v - l) / std::min(l,1.0-l);
}

inline color_hsl::color_hsl(color_lrgb z)
{
    num_t xmax = std::max(std::max(z.r,z.g),z.b);
    num_t xmin = std::min(std::min(z.r,z.g),z.b);
    num_t c = xmax - xmin;
    l = (xmax + xmin) / 2.0;
    s = (l == 0.0 || l == 1.0) ? 0.0 : (xmax-l)/std::min(l,1.0-l);
    h = _hue_from_rgb(c,xmax,z.r,z.g,z.b);
}

inline color_hsv::color_hsv(color_hsl z)
{
    h = z.h;
    v = z.l + z.s*std::min(z.l,1.0-z.l);
    s = (v == 0.0) ? 0.0 : 2.0*(1.0 - z.l/v);
}

inline color_hsv::color_hsv(color_lrgb z)
{
    num_t xmax = std::max(std::max(z.r,z.g),z.b);
    num_t xmin = std::min(std::min(z.r,z.g),z.b);
    num_t c = xmax - xmin;
    v = xmax;
    s = (v == 0) ? 0.0 : c/v;
    h = _hue_from_rgb(c,v,z.r,z.g,z.b);
}

inline color_xyz::color_xyz(color_lrgb w)
{
    x = 0.4124*w.r + 0.3576*w.g + 0.1805*w.b;
    y = 0.2126*w.r + 0.7152*w.g + 0.0722*w.b;
    z = 0.0193*w.r + 0.1192*w.g + 0.9505*w.b;
}

inline color_xyz::color_xyz(color_oklab w)
{
    num_t l = 0.9999999984505196*w.L + 0.39633779217376774*w.a + 0.2158037580607588*w.b;
    num_t m = 1.0000000088817607*w.L - 0.10556134232365633*w.a - 0.0638541747717059*w.b;
    num_t s = 1.0000000546724108*w.L - 0.08948418209496574*w.a - 1.2914855378640917*w.b;
    l = l*l*l;
    m = m*m*m;
    s = s*s*s;
    x = 1.2270138511035211*l - 0.5577999806518222*m + 0.2812561489664678*s;
    y = -0.0405801784232806*l + 1.11225686961683*m - 0.07167667866560119*s;
    z = -0.0763812845057069*l - 0.4214819784180127*m + 1.5861632204407947*s;
}

inline color_oklab::color_oklab(color_xyz w)
{
    num_t l = 0.8189330101*w.x + 0.3618667424*w.y - 0.1288597137*w.z;
    num_t m = 0.0329845436*w.x + 0.9293118715*w.y + 0.0361456387*w.z;
    num_t s = 0.0482003018*w.x + 0.2643662691*w.y + 0.6338517070*w.z;
    l = std::cbrt(l);
    m = std::cbrt(m);
    s = std::cbrt(s);
    L = 0.2104542553*l + 0.7936177850*m - 0.0040720468*s;
    a = 1.9779984951*l - 2.4285922050*m + 0.4505937099*s;
    b = 0.0259040371*l + 0.7827717662*m - 0.8086757660*s;
}

inline color_oklab::color_oklab(color_lrgb z)
{
    num_t l = 0.4122214708*z.r + 0.5363325363*z.g + 0.0514459929*z.b;
    num_t m = 0.2119034982*z.r + 0.6806995451*z.g + 0.1073969566*z.b;
    num_t s = 0.0883024619*z.r + 0.2817188376*z.g + 0.6299787005*z.b;
    l = std::cbrt(l);
    m = std::cbrt(m);
    s = std::cbrt(s);
    L = 0.2104542553*l + 0.2104542553*m - 0.0040720468*s;
    a = 1.9779984951*l - 2.4285922050*m + 0.4505937099*s;
    b = 0.0259040371*l + 0.7827717662*m - 0.8086757660*s;
}

}
