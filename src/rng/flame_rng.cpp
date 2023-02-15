#include "flame_rng.hpp"

namespace tkoz::flame
{

template class FlameRNG<float,u32,4>;
template class FlameRNG<float,u32,6>;
template class FlameRNG<float,u32,8>;

template class FlameRNG<double,u32,4>;
template class FlameRNG<double,u32,6>;
template class FlameRNG<double,u32,8>;

template class FlameRNG<float,u64,4>;
template class FlameRNG<float,u64,6>;
template class FlameRNG<float,u64,8>;

template class FlameRNG<double,u64,4>;
template class FlameRNG<double,u64,6>;
template class FlameRNG<double,u64,8>;

}
