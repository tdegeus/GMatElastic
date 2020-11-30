/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_ARRAY_HPP
#define GMATELASTIC_CARTESIAN3D_ARRAY_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

template <size_t N>
inline Array<N>::Array(const std::array<size_t, N>& shape)
{
    this->init(shape);
    m_allSet = false;
    m_type = xt::ones<size_t>(m_shape) * Type::Unset;
    m_index = xt::empty<size_t>(m_shape);
}

template <size_t N>
inline Array<N>::Array(const std::array<size_t, N>& shape, double K, double G)
{
    this->init(shape);
    m_allSet = false;
    m_type = xt::ones<size_t>(m_shape) * Type::Elastic;
    m_index = xt::arange<size_t>(m_size).reshape(m_shape);

    for (size_t i = 0; i < m_size; ++i) {
        m_Elastic.push_back(Elastic(K, G));
    }
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::K() const
{
    GMATELASTIC_ASSERT(m_allSet);
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            ret.data()[i] = m_Elastic[m_index.data()[i]].K();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::G() const
{
    GMATELASTIC_ASSERT(m_allSet);
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            ret.data()[i] = m_Elastic[m_index.data()[i]].G();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::type() const
{
    GMATELASTIC_ASSERT(m_allSet);
    return m_type;
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::isElastic() const
{
    GMATELASTIC_ASSERT(m_allSet);
    xt::xtensor<size_t, N> ret = xt::where(xt::equal(m_type, Type::Elastic), 1ul, 0ul);
    return ret;
}

template <size_t N>
inline void Array<N>::check() const
{
    if (xt::any(xt::equal(m_type, Type::Unset))) {
        throw std::runtime_error("Points without material found");
    }
}

template <size_t N>
inline void Array<N>::checkAllSet()
{
    if (xt::any(xt::equal(m_type, Type::Unset))) {
        m_allSet = false;
    }
    else {
        m_allSet = true;
    }
}

template <size_t N>
inline void Array<N>::setElastic(const xt::xtensor<size_t, N>& I, double K, double G)
{
    GMATELASTIC_ASSERT(m_type.shape() == I.shape());
    GMATELASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATELASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    for (size_t i = 0; i < m_size; ++i) {
        if (I.data()[i] == 1ul) {
            m_type.data()[i] = Type::Elastic;
            m_index.data()[i] = m_Elastic.size();
            m_Elastic.push_back(Elastic(K, G));
        }
    }

    this->checkAllSet();
}

template <size_t N>
inline void Array<N>::setStrain(const xt::xtensor<double, N + 2>& A)
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].setStrainIterator(&A.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::strain(xt::xtensor<double, N + 2>& A) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].strainIterator(&A.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::stress(xt::xtensor<double, N + 2>& A) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].stressIterator(&A.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::tangent(xt::xtensor<double, N + 4>& A) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor4));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        auto c = xt::adapt(&A.data()[i * m_stride_tensor4], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].tangent(c);
            break;
        }
    }
}

template <size_t N>
inline xt::xtensor<double, N + 2> Array<N>::Strain() const
{
    xt::xtensor<double, N + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->strain(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N + 2> Array<N>::Stress() const
{
    xt::xtensor<double, N + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->stress(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N + 4> Array<N>::Tangent() const
{
    xt::xtensor<double, N + 4> ret = xt::empty<double>(m_shape_tensor4);
    this->tangent(ret);
    return ret;
}

template <size_t N>
inline auto Array<N>::getElastic(const std::array<size_t, N>& index) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(m_type[index] == Type::Elastic);
    return m_Elastic[m_index[index]];
}

template <size_t N>
inline auto* Array<N>::refElastic(const std::array<size_t, N>& index)
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(m_type[index] == Type::Elastic);
    return &m_Elastic[m_index[index]];
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
