/*
Representation of the variation information, the function and parameters parser
*/

#pragma once

#include <cstdlib>

namespace tkoz
{
namespace flame
{

template <typename num_t, size_t dims>
using VarFunc = std::function<void(IterState<num_t,dims>&,const num_t*)>;

template <typename num_t, size_t dims>
using VarParse = std::function<void(const Json&,std::vector<num_t>&)>;

// representation of variations data
template <typename num_t, size_t dims>
class VarInfo
{
private:
    // function pointer for calculating the variation
    VarFunc<num_t,dims> func;
    // params parser, use nullptr for default of weight only
    VarParse<num_t,dims> params;
public:
    // default constructable for use with unordered_map operator[]
    VarInfo(VarFunc<num_t,dims> func = nullptr,
            VarParse<num_t,dims> params = nullptr):
        func(func),params(params) {}
    inline VarFunc<num_t,dims> getFPtr() const
    {
        return func;
    }
    inline VarParse<num_t,dims> getPPtr() const
    {
        return params;
    }
};

template <typename num_t, size_t dims>
using VarData = std::unordered_map<std::string,VarInfo<num_t,dims>>;

}
}
