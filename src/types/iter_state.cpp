#include "iter_state.hpp"

namespace tkoz::flame
{
#if INSTANTIATE_TEMPLATES
template class IterState<float,1>;
template class IterState<float,2>;
template class IterState<float,3>;
template class IterState<float,4>;
template class IterState<float,5>;

template class IterState<double,1>;
template class IterState<double,2>;
template class IterState<double,3>;
template class IterState<double,4>;
template class IterState<double,5>;
#endif
}
