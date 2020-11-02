/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_ARRAY_HPP
#define GMATELASTIC_CARTESIAN3D_ARRAY_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

template <size_t rank>
inline Array<rank>::Array(const std::array<size_t, rank>& shape) : m_shape(shape)
{
    m_allSet = false;
    size_t nd = m_ndim;
    std::copy(shape.begin(), shape.end(), m_shape_tensor2.begin());
    std::copy(shape.begin(), shape.end(), m_shape_tensor4.begin());
    std::fill(m_shape_tensor2.begin() + rank, m_shape_tensor2.end(), nd);
    std::fill(m_shape_tensor4.begin() + rank, m_shape_tensor4.end(), nd);
    m_size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    m_type = xt::ones<size_t>(m_shape) * Type::Unset;
    m_index = xt::empty<size_t>(m_shape);
}

template <size_t rank>
inline Array<rank>::Array(const std::array<size_t, rank>& shape, double K, double G) : m_shape(shape)
{
    m_allSet = false;
    size_t nd = m_ndim;
    std::copy(shape.begin(), shape.end(), m_shape_tensor2.begin());
    std::copy(shape.begin(), shape.end(), m_shape_tensor4.begin());
    std::fill(m_shape_tensor2.begin() + rank, m_shape_tensor2.end(), nd);
    std::fill(m_shape_tensor4.begin() + rank, m_shape_tensor4.end(), nd);
    m_size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    m_type = xt::ones<size_t>(m_shape) * Type::Elastic;
    m_index = xt::arange<size_t>(m_size).reshape(m_shape);

    for (size_t i = 0; i < m_size; ++i) {
        m_Elastic.push_back(Elastic(K, G));
    }
}

template <size_t rank>
inline std::array<size_t, rank> Array<rank>::shape() const
{
    return m_shape;
}

template <size_t rank>
inline xt::xtensor<double, rank> Array<rank>::K() const
{
    GMATELASTIC_ASSERT(m_allSet);
    xt::xtensor<double, rank> ret = xt::empty<double>(m_shape);

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

template <size_t rank>
inline xt::xtensor<double, rank> Array<rank>::G() const
{
    GMATELASTIC_ASSERT(m_allSet);
    xt::xtensor<double, rank> ret = xt::empty<double>(m_shape);

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

template <size_t rank>
inline xt::xtensor<double, rank + 2> Array<rank>::I2() const
{
    xt::xtensor<double, rank + 2> ret = xt::empty<double>(m_shape_tensor2);

    #pragma omp parallel
    {
        Tensor2 unit = Cartesian3d::I2();
        size_t stride = m_ndim * m_ndim;

        #pragma omp for
        for (size_t i = 0; i < m_size; ++i) {
            auto view = xt::adapt(&ret.data()[i * stride], xt::xshape<m_ndim, m_ndim>());
            xt::noalias(view) = unit;
        }
    }

    return ret;
}

template <size_t rank>
inline xt::xtensor<double, rank + 4> Array<rank>::II() const
{
    xt::xtensor<double, rank + 4> ret = xt::empty<double>(m_shape_tensor4);

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::II();
        size_t stride = m_ndim * m_ndim * m_ndim * m_ndim;

        #pragma omp for
        for (size_t i = 0; i < m_size; ++i) {
            auto view = xt::adapt(&ret.data()[i * stride], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
            xt::noalias(view) = unit;
        }
    }

    return ret;
}

template <size_t rank>
inline xt::xtensor<double, rank + 4> Array<rank>::I4() const
{
    xt::xtensor<double, rank + 4> ret = xt::empty<double>(m_shape_tensor4);

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4();
        size_t stride = m_ndim * m_ndim * m_ndim * m_ndim;

        #pragma omp for
        for (size_t i = 0; i < m_size; ++i) {
            auto view = xt::adapt(&ret.data()[i * stride], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
            xt::noalias(view) = unit;
        }
    }

    return ret;
}

template <size_t rank>
inline xt::xtensor<double, rank + 4> Array<rank>::I4rt() const
{
    xt::xtensor<double, rank + 4> ret = xt::empty<double>(m_shape_tensor4);

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4rt();
        size_t stride = m_ndim * m_ndim * m_ndim * m_ndim;

        #pragma omp for
        for (size_t i = 0; i < m_size; ++i) {
            auto view = xt::adapt(&ret.data()[i * stride], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
            xt::noalias(view) = unit;
        }
    }

    return ret;
}

template <size_t rank>
inline xt::xtensor<double, rank + 4> Array<rank>::I4s() const
{
    xt::xtensor<double, rank + 4> ret = xt::empty<double>(m_shape_tensor4);

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4s();
        size_t stride = m_ndim * m_ndim * m_ndim * m_ndim;

        #pragma omp for
        for (size_t i = 0; i < m_size; ++i) {
            auto view = xt::adapt(&ret.data()[i * stride], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
            xt::noalias(view) = unit;
        }
    }

    return ret;
}

template <size_t rank>
inline xt::xtensor<double, rank + 4> Array<rank>::I4d() const
{
    xt::xtensor<double, rank + 4> ret = xt::empty<double>(m_shape_tensor4);

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4d();
        size_t stride = m_ndim * m_ndim * m_ndim * m_ndim;

        #pragma omp for
        for (size_t i = 0; i < m_size; ++i) {
            auto view = xt::adapt(&ret.data()[i * stride], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
            xt::noalias(view) = unit;
        }
    }

    return ret;
}

template <size_t rank>
inline xt::xtensor<size_t, rank> Array<rank>::type() const
{
    GMATELASTIC_ASSERT(m_allSet);
    return m_type;
}

template <size_t rank>
inline xt::xtensor<size_t, rank> Array<rank>::isElastic() const
{
    GMATELASTIC_ASSERT(m_allSet);
    xt::xtensor<size_t, rank> ret = xt::where(xt::equal(m_type, Type::Elastic), 1ul, 0ul);
    return ret;
}

template <size_t rank>
inline void Array<rank>::check() const
{
    if (xt::any(xt::equal(m_type, Type::Unset))) {
        throw std::runtime_error("Points without material found");
    }
}

template <size_t rank>
inline void Array<rank>::checkAllSet()
{
    if (xt::any(xt::equal(m_type, Type::Unset))) {
        m_allSet = false;
    }
    else {
        m_allSet = true;
    }
}

template <size_t rank>
inline void Array<rank>::setElastic(const xt::xtensor<size_t, rank>& I, double K, double G)
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

template <size_t rank>
inline void Array<rank>::setStrain(const xt::xtensor<double, rank + 2>& A)
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].setStrainIterator(&A.data()[i * m_ndim * m_ndim]);
            break;
        }
    }
}

template <size_t rank>
inline void Array<rank>::stress(xt::xtensor<double, rank + 2>& A) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].stressIterator(&A.data()[i * m_ndim * m_ndim]);
            break;
        }
    }
}

template <size_t rank>
inline void Array<rank>::tangent(xt::xtensor<double, rank + 4>& A) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor4));
    size_t stride = m_ndim * m_ndim * m_ndim * m_ndim;

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        auto c = xt::adapt(&A.data()[i * stride], xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());
        switch (m_type.data()[i]) {
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].tangent(c);
            break;
        }
    }
}

template <size_t rank>
inline xt::xtensor<double, rank + 2> Array<rank>::Stress() const
{
    xt::xtensor<double, rank + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->stress(ret);
    return ret;
}

template <size_t rank>
inline xt::xtensor<double, rank + 4> Array<rank>::Tangent() const
{
    xt::xtensor<double, rank + 4> ret = xt::empty<double>(m_shape_tensor4);
    this->tangent(ret);
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
