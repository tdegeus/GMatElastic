/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTIC_CARTESIAN3D_HPP
#define GMATELASTIC_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

template <class T, class U>
inline void epseq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(2.0 / 3.0);
}

template <class T>
inline auto Epseq(const T& A) ->
    typename GMatTensor::detail::allocate<xt::get_rank<T>::value - 2, T>::type
{
    return xt::eval(std::sqrt(2.0 / 3.0) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

template <class T, class U>
inline void sigeq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(1.5);
}

template <class T>
inline auto Sigeq(const T& A) ->
    typename GMatTensor::detail::allocate<xt::get_rank<T>::value - 2, T>::type
{
    return xt::eval(std::sqrt(1.5) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
