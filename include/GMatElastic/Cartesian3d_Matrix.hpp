/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_MATRIX_HPP
#define GMATELASTIC_CARTESIAN3D_MATRIX_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

inline Matrix::Matrix(size_t nelem, size_t nip) : m_nelem(nelem), m_nip(nip)
{
    m_type = xt::ones<size_t>({m_nelem, m_nip}) * Type::Unset;
    m_index = xt::empty<size_t>({m_nelem, m_nip});
    m_allSet = false;
}

inline Matrix::Matrix(size_t nelem, size_t nip, double K, double G) : m_nelem(nelem), m_nip(nip)
{
    m_type = xt::ones<size_t>({m_nelem, m_nip}) * Type::Elastic;
    m_index = xt::ones<size_t>({m_nelem, m_nip}) * m_Elastic.size();
    m_allSet = true;

    m_Elastic.push_back(Elastic(K, G));
}

inline size_t Matrix::ndim() const
{
    return m_ndim;
}

inline size_t Matrix::nelem() const
{
    return m_nelem;
}

inline size_t Matrix::nip() const
{
    return m_nip;
}

inline xt::xtensor<size_t,2> Matrix::type() const
{
    return m_type;
}

inline xt::xtensor<double,2> Matrix::K() const
{
    GMATELASTIC_ASSERT(m_allSet);

    xt::xtensor<double,2> out = xt::empty<double>({m_nelem, m_nip});

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t q = 0; q < m_nip; ++q) {

            switch (m_type(e, q)) {
            case Type::Elastic:
                out(e, q) = m_Elastic[m_index(e, q)].K();
                break;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,2> Matrix::G() const
{
    GMATELASTIC_ASSERT(m_allSet);

    xt::xtensor<double,2> out = xt::empty<double>({m_nelem, m_nip});

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t q = 0; q < m_nip; ++q) {

            switch (m_type(e, q)) {
            case Type::Elastic:
                out(e, q) = m_Elastic[m_index(e, q)].G();
                break;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,4> Matrix::I2() const
{
    xt::xtensor<double,4> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

    #pragma omp parallel
    {
        Tensor2 unit = Cartesian3d::I2();

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t q = 0; q < m_nip; ++q) {
                auto view = xt::adapt(&out(e, q, 0, 0), xt::xshape<m_ndim, m_ndim>());
                xt::noalias(view) = unit;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,6> Matrix::II() const
{
    xt::xtensor<double,6> out =
        xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::II();

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t q = 0; q < m_nip; ++q) {

                auto view =
                    xt::adapt(&out(e, q, 0, 0, 0, 0), xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());

                xt::noalias(view) = unit;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,6> Matrix::I4() const
{
    xt::xtensor<double,6> out =
        xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4();

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t q = 0; q < m_nip; ++q) {

                auto view =
                    xt::adapt(&out(e, q, 0, 0, 0, 0), xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());

                xt::noalias(view) = unit;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,6> Matrix::I4rt() const
{
    xt::xtensor<double,6> out =
        xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4rt();

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t q = 0; q < m_nip; ++q) {

                auto view =
                    xt::adapt(&out(e, q, 0, 0, 0, 0), xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());

                xt::noalias(view) = unit;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,6> Matrix::I4s() const
{
    xt::xtensor<double,6> out =
        xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4s();

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t q = 0; q < m_nip; ++q) {

                auto view =
                    xt::adapt(&out(e, q, 0, 0, 0, 0), xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());

                xt::noalias(view) = unit;
            }
        }
    }

    return out;
}

inline xt::xtensor<double,6> Matrix::I4d() const
{
    xt::xtensor<double,6> out =
        xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

    #pragma omp parallel
    {
        Tensor4 unit = Cartesian3d::I4d();

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t q = 0; q < m_nip; ++q) {

                auto view =
                    xt::adapt(&out(e, q, 0, 0, 0, 0), xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());

                xt::noalias(view) = unit;
            }
        }
    }

    return out;
}

inline void Matrix::check() const
{
    if (xt::any(xt::equal(m_type, Type::Unset))) {
        throw std::runtime_error("Points without material found");
    }
}

inline void Matrix::checkAllSet()
{
    if (xt::any(xt::equal(m_type, Type::Unset))) {
        m_allSet = false;
    }
    else {
        m_allSet = true;
    }
}

inline void Matrix::setElastic(const xt::xtensor<size_t,2>& I, double K, double G)
{
    GMATELASTIC_ASSERT(m_type.shape() == I.shape());
    GMATELASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATELASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    m_type = xt::where(xt::equal(I, 1ul), Type::Elastic, m_type);
    m_index = xt::where(xt::equal(I, 1ul), m_Elastic.size(), m_index);
    this->checkAllSet();
    m_Elastic.push_back(Elastic(K, G));
}

inline void Matrix::stress(const xt::xtensor<double,4>& a_Eps, xt::xtensor<double,4>& a_Sig) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(
        a_Eps.shape() ==
        std::decay_t<decltype(a_Eps)>::shape_type({m_nelem, m_nip, m_ndim, m_ndim}));
    GMATELASTIC_ASSERT(a_Eps.shape() == a_Sig.shape());

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t q = 0; q < m_nip; ++q) {

            auto Eps = xt::adapt(&a_Eps(e, q, 0, 0), xt::xshape<m_ndim, m_ndim>());
            auto Sig = xt::adapt(&a_Sig(e, q, 0, 0), xt::xshape<m_ndim, m_ndim>());

            switch (m_type(e, q)) {
            case Type::Elastic:
                m_Elastic[m_index(e, q)].stress(Eps, Sig);
                break;
            }
        }
    }
}

inline void Matrix::tangent(
    const xt::xtensor<double,4>& a_Eps,
          xt::xtensor<double,4>& a_Sig,
          xt::xtensor<double,6>& a_C) const
{
    GMATELASTIC_ASSERT(m_allSet);
    GMATELASTIC_ASSERT(
        a_Eps.shape() ==
        std::decay_t<decltype(a_Eps)>::shape_type({m_nelem, m_nip, m_ndim, m_ndim}));
    GMATELASTIC_ASSERT(a_Eps.shape() == a_Sig.shape());
    GMATELASTIC_ASSERT(
        a_C.shape() ==
        std::decay_t<decltype(a_C)>::shape_type({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim}));

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t q = 0; q < m_nip; ++q) {

            auto Eps = xt::adapt(&a_Eps(e, q, 0, 0), xt::xshape<m_ndim, m_ndim>());
            auto Sig = xt::adapt(&a_Sig(e, q, 0, 0), xt::xshape<m_ndim, m_ndim>());

            auto C =
                xt::adapt(&a_C(e, q, 0, 0, 0, 0), xt::xshape<m_ndim, m_ndim, m_ndim, m_ndim>());

            switch (m_type(e, q)) {
            case Type::Elastic:
                m_Elastic[m_index(e, q)].tangent(Eps, Sig, C);
                break;
            }
        }
    }
}

inline xt::xtensor<double,4> Matrix::Stress(const xt::xtensor<double,4>& Eps) const
{
    xt::xtensor<double,4> Sig = xt::empty<double>(Eps.shape());
    this->stress(Eps, Sig);
    return Sig;
}

inline std::tuple<xt::xtensor<double,4>, xt::xtensor<double,6>>
Matrix::Tangent(const xt::xtensor<double,4>& Eps) const
{
    xt::xtensor<double,4> Sig = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});
    xt::xtensor<double,6> C = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});
    this->tangent(Eps, Sig, C);
    return std::make_tuple(Sig, C);
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
