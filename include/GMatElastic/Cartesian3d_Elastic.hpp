/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTIC_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

inline double Elastic::K() const
{
    return m_K;
}

inline double Elastic::G() const
{
    return m_G;
}

template <class T>
inline void Elastic::setStrain(const T& a)
{
    GMATELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->setStrainIterator(a.cbegin());
}

template <class T>
inline void Elastic::setStrainIterator(const T& begin)
{
    std::copy(begin, begin + 9, m_Eps.begin());

    double epsm = detail::trace(m_Eps) / 3.0;

    m_Sig[0] = (m_K - m_G) * epsm + m_G * m_Eps[0];
    m_Sig[1] = m_G * m_Eps[1];
    m_Sig[2] = m_G * m_Eps[2];
    m_Sig[3] = m_G * m_Eps[3];
    m_Sig[4] = (m_K - m_G) * epsm + m_G * m_Eps[4];
    m_Sig[5] = m_G * m_Eps[5];
    m_Sig[6] = m_G * m_Eps[6];
    m_Sig[7] = m_G * m_Eps[7];
    m_Sig[8] = (m_K - m_G) * epsm + m_G * m_Eps[8];
}

template <class T>
inline void Elastic::stress(T& a) const
{
    GMATELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->stressIterator(a.begin());
}

template <class T>
inline void Elastic::stressIterator(const T& begin) const
{
    std::copy(m_Sig.begin(), m_Sig.end(), begin);
}

inline Tensor2 Elastic::Stress() const
{
    auto ret = Tensor2::from_shape({3, 3});
    this->stressIterator(ret.begin());
    return ret;
}

template <class T>
inline void Elastic::tangent(T& C) const
{
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    xt::noalias(C) = m_K * II + 2.0 * m_G * I4d;
}

inline Tensor4 Elastic::Tangent() const
{
    auto ret = Tensor4::from_shape({3, 3, 3, 3});
    this->tangent(ret);
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
