#include "renderer.hpp"

namespace tkoz::flame
{
#if INSTANTIATE_TEMPLATES
template class RendererBasic<float,u32>;
template class RendererBasic<double,u32>;
template class RendererBasic<float,u64>;
template class RendererBasic<double,u64>;
#endif
}
