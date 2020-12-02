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
inline void Elastic::setStrainIterator(const T* arg)
{
    using namespace GMatTensor::Cartesian3d::pointer;
    std::copy(arg, arg + 9, m_Eps.begin());

    double epsm = trace(&m_Eps[0]) / 3.0;

    m_Sig[0] = (3.0 * m_K - 2.0 * m_G) * epsm + 2.0 * m_G * m_Eps[0];
    m_Sig[1] = 2.0 * m_G * m_Eps[1];
    m_Sig[2] = 2.0 * m_G * m_Eps[2];
    m_Sig[3] = 2.0 * m_G * m_Eps[3];
    m_Sig[4] = (3.0 * m_K - 2.0 * m_G) * epsm + 2.0 * m_G * m_Eps[4];
    m_Sig[5] = 2.0 * m_G * m_Eps[5];
    m_Sig[6] = 2.0 * m_G * m_Eps[6];
    m_Sig[7] = 2.0 * m_G * m_Eps[7];
    m_Sig[8] = (3.0 * m_K - 2.0 * m_G) * epsm + 2.0 * m_G * m_Eps[8];
}

template <class T>
inline void Elastic::strainIterator(T* ret) const
{
    std::copy(m_Eps.begin(), m_Eps.end(), ret);
}

template <class T>
inline void Elastic::stressIterator(T* ret) const
{
    std::copy(m_Sig.begin(), m_Sig.end(), ret);
}

template <class T>
inline void Elastic::tangentIterator(T* ret) const
{
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    auto C = m_K * II + 2.0 * m_G * I4d;
    std::copy(C.cbegin(), C.cend(), ret);
}

template <class T>
inline void Elastic::setStrain(const T& a)
{
    GMATELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->setStrainIterator(a.data());
}

template <class T>
inline void Elastic::strain(T& a) const
{
    GMATELASTOPLASTICQPOT_ASSERT(xt::has_shape(a, {3, 3}));
    return this->strainIterator(a.data());
}

template <class T>
inline void Elastic::stress(T& a) const
{
    GMATELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->stressIterator(a.data());
}

template <class T>
inline void Elastic::tangent(T& a) const
{
    GMATELASTIC_ASSERT(xt::has_shape(a, {3, 3, 3, 3}));
    return this->stressIterator(a.data());
}

inline xt::xtensor<double, 2> Elastic::Strain() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->strainIterator(ret.data());
    return ret;
}

inline xt::xtensor<double, 2> Elastic::Stress() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->stressIterator(ret.data());
    return ret;
}

inline xt::xtensor<double, 4> Elastic::Tangent() const
{
    xt::xtensor<double, 4> ret = xt::empty<double>({3, 3, 3, 3});
    this->tangentIterator(ret.data());
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
