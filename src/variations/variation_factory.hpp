#pragma once

#include <cstdlib>

#include "../utils/json_small.hpp"
#include "variation_base.hpp"

// macros for using SFINAE
#define ENABLE_IF2(COND) template <typename RET> \
    typename std::enable_if<(COND),RET>::type

namespace tkoz::flame::vars
{

// do not include 2d variations with below 2 dimensions
template <typename num_t, size_t dims>
ENABLE_IF2(dims<2) Variation<num_t,dims>::parseVariation2d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put 2d variation factory here
template <typename num_t, size_t dims>
ENABLE_IF2(dims>=2) Variation<num_t,dims>::parseVariation2d(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "swirl")
        return new Swirl<num_t,dims>(json);
    if (name == "horseshoe")
        return new Horseshoe<num_t,dims>(json);
    if (name == "polar")
        return new Polar<num_t,dims>(json);
    if (name == "polar2")
        return new Polar2<num_t,dims>(json);
    if (name == "handkerchief")
        return new Handkerchief<num_t,dims>(json);
    if (name == "heart")
        return new Heart<num_t,dims>(json);
    if (name == "disc")
        return new Disc<num_t,dims>(json);
    if (name == "disc2")
        return new Disc2<num_t,dims>(json);
    if (name == "waves")
        return new Waves<num_t,dims>(json);
    if (name == "fan")
        return new Fan<num_t,dims>(json);
    if (name == "rings")
        return new Rings<num_t,dims>(json);
    if (name == "spiral")
        return new Spiral<num_t,dims>(json);
    if (name == "hyperbolic")
        return new Hyperbolic<num_t,dims>(json);
    if (name == "diamond")
        return new Diamond<num_t,dims>(json);
    if (name == "ex")
        return new Ex<num_t,dims>(json);
    if (name == "julia")
        return new Julia<num_t,dims>(json);
    if (name == "exponential")
        return new Exponential<num_t,dims>(json);
    if (name == "power")
        return new Power<num_t,dims>(json);
    if (name == "cosine")
        return new Cosine<num_t,dims>(json);
    if (name == "blob")
        return new Blob<num_t,dims>(json);
    if (name == "pdj")
        return new PDJ<num_t,dims>(json);
    if (name == "cylinder")
        return new Cylinder<num_t,dims>(json);
    if (name == "perspective")
        return new Perspective<num_t,dims>(json);
    if (name == "julian")
        return new JuliaN<num_t,dims>(json);
    if (name == "juliascope")
        return new JuliaScope<num_t,dims>(json);
    if (name == "radial_blur")
        return new RadialBlur<num_t,dims>(json);
    if (name == "pie")
        return new Pie<num_t,dims>(json);
    if (name == "ngon")
        return new NGon<num_t,dims>(json);
    if (name == "curl")
        return new Curl<num_t,dims>(json);
    if (name == "arch")
        return new Arch<num_t,dims>(json);
    if (name == "tangent")
        return new Tangent<num_t,dims>(json);
    if (name == "rays")
        return new Rays<num_t,dims>(json);
    if (name == "blade")
        return new Blade<num_t,dims>(json);
    if (name == "secant")
        return new Secant<num_t,dims>(json);
    if (name == "twintrian")
        return new Twintrian<num_t,dims>(json);
    if (name == "cross")
        return new Cross<num_t,dims>(json);
    if (name == "exp")
        return new Exp<num_t,dims>(json);
    if (name == "log")
        return new Log<num_t,dims>(json);
    if (name == "sin")
        return new Sin<num_t,dims>(json);
    if (name == "cos")
        return new Cos<num_t,dims>(json);
    if (name == "tan")
        return new Tan<num_t,dims>(json);
    if (name == "sec")
        return new Sec<num_t,dims>(json);
    if (name == "csc")
        return new Csc<num_t,dims>(json);
    if (name == "cot")
        return new Cot<num_t,dims>(json);
    if (name == "sinh")
        return new Sinh<num_t,dims>(json);
    if (name == "cosh")
        return new Cosh<num_t,dims>(json);
    if (name == "tanh")
        return new Tanh<num_t,dims>(json);
    if (name == "sech")
        return new Sech<num_t,dims>(json);
    if (name == "csch")
        return new Csch<num_t,dims>(json);
    if (name == "coth")
        return new Coth<num_t,dims>(json);
    if (name == "auger")
        return new Auger<num_t,dims>(json);
    if (name == "flux")
        return new Flux<num_t,dims>(json);
    if (name == "mobius")
        return new Mobius<num_t,dims>(json);
    if (name == "scry")
        return new Scry<num_t,dims>(json);
    if (name == "split")
        return new Split<num_t,dims>(json);
    if (name == "stripes")
        return new Stripes<num_t,dims>(json);
    if (name == "wedge")
        return new Wedge<num_t,dims>(json);
    if (name == "wedge_julia")
        return new WedgeJulia<num_t,dims>(json);
    if (name == "wedge_sph")
        return new WedgeSph<num_t,dims>(json);
    if (name == "whorl")
        return new Whorl<num_t,dims>(json);
    if (name == "supershape")
        return new Supershape<num_t,dims>(json);
    if (name == "flower")
        return new Flower<num_t,dims>(json);
    if (name == "conic")
        return new Conic<num_t,dims>(json);
    if (name == "parabola")
        return new Parabola<num_t,dims>(json);
    if (name == "bipolar")
        return new Bipolar<num_t,dims>(json);
    if (name == "boarders")
        return new Boarders<num_t,dims>(json);
    if (name == "butterfly")
        return new Butterfly<num_t,dims>(json);
    if (name == "cell")
        return new Cell<num_t,dims>(json);
    if (name == "cpow")
        return new CPow<num_t,dims>(json);
    if (name == "curve")
        return new Curve<num_t,dims>(json);
    if (name == "edisc")
        return new EDisc<num_t,dims>(json);
    if (name == "elliptic")
        return new Elliptic<num_t,dims>(json);
    if (name == "escher")
        return new Escher<num_t,dims>(json);
    if (name == "foci")
        return new Foci<num_t,dims>(json);
    if (name == "lazysusan")
        return new LazySusan<num_t,dims>(json);
    if (name == "loonie")
        return new Loonie<num_t,dims>(json);
    if (name == "oscope")
        return new OScope<num_t,dims>(json);
    if (name == "popcorn")
        return new Popcorn<num_t,dims>(json);
    return nullptr;
}

// do not include 3d variations with below 3 dimensions
template <typename num_t, size_t dims>
ENABLE_IF2(dims<3) Variation<num_t,dims>::parseVariation3d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put 3d variation factory here
template <typename num_t, size_t dims>
ENABLE_IF2(dims>=3) Variation<num_t,dims>::parseVariation3d(const Json& json)
{
    (void)json;
    return nullptr;
}

// put nd variation factory here
template <typename num_t, size_t dims>
Variation<num_t,dims> *Variation<num_t,dims>::parseVariationNd(const Json& json)
{
    std::string name = json["name"].stringValue();
    if (name == "linear")
        return new Linear<num_t,dims>(json);
    if (name == "sinusoidal")
        return new Sinusoidal<num_t,dims>(json);
    if (name == "spherical")
        return new Spherical<num_t,dims>(json);
    if (name == "bent")
        return new Bent<num_t,dims>(json);
    if (name == "rectangles")
        return new Rectangles<num_t,dims>(json);
    if (name == "fisheye")
        return new Fisheye<num_t,dims>(json);
    if (name == "bubble")
        return new Bubble<num_t,dims>(json);
    if (name == "noise")
        return new Noise<num_t,dims>(json);
    if (name == "blur")
        return new Blur<num_t,dims>(json);
    if (name == "gaussian_blur")
        return new GaussianBlur<num_t,dims>(json);
    if (name == "square_noise")
        return new SquareNoise<num_t,dims>(json);
    if (name == "separation")
        return new Separation<num_t,dims>(json);
    if (name == "splits")
        return new Splits<num_t,dims>(json);
    if (name == "pre_blur")
        return new PreBlur<num_t,dims>(json);
    if (name == "modulus")
        return new Modulus<num_t,dims>(json);
    if (name == "spherical_p")
        return new SphericalP<num_t,dims>(json);
    if (name == "unit_sphere")
        return new UnitSphere<num_t,dims>(json);
    if (name == "unit_sphere_p")
        return new UnitSphereP<num_t,dims>(json);
    if (name == "unit_cube")
        return new UnitCube<num_t,dims>(json);
    return nullptr;
}

template <typename num_t, size_t dims>
Variation<num_t,dims> *Variation<num_t,dims>::parseVariation(const Json& json)
{
    std::string name = json["name"].stringValue();
    Variation<num_t,dims> *ret = nullptr;
    ret = parseVariation2d(json);
    if (ret) return ret;
    ret = parseVariation3d(json);
    if (ret) return ret;
    ret = parseVariationNd(json);
    if (ret) return ret;
    throw std::runtime_error("unknown variation: "+name);
}

}

#undef ENABLE_IF2
