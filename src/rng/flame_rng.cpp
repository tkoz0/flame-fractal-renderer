#include "flame_rng.hpp"

namespace tkoz::flame
{

// static thread local rng object requires initialization
template <> thread_local
Isaac<hist_t,4> FlameRNG<num_t,hist_t,4>::state = Isaac<hist_t,4>();

} // namespace tkoz::flame
