#include "affine.hpp"

namespace tkoz::flame
{
#if INSTANTIATE_TEMPLATES
template class Affine<float,1>;
template class Affine<float,2>;
template class Affine<float,3>;
template class Affine<float,4>;
template class Affine<float,5>;

template class Affine<double,1>;
template class Affine<double,2>;
template class Affine<double,3>;
template class Affine<double,4>;
template class Affine<double,5>;
#endif
}
