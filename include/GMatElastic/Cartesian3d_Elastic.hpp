/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

================================================================================================= */

#ifndef GMATELASTIC_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTIC_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::K() const
{
  return m_K;
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::G() const
{
  return m_G;
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline void Elastic::stress(const T2& Eps, T&& Sig) const
{
  auto I     = Cartesian3d::I();
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps / 3.0 * I;
  xt::noalias(Sig) = m_K * treps * I + 2.0 * m_G * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline T2 Elastic::Stress(const T2& Eps) const
{
  T2 Sig;
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

template <class T, class S>
inline void Elastic::tangent(const T2& Eps, T&& Sig, S&& C) const
{
  auto I     = Cartesian3d::I();
  auto II    = Cartesian3d::II();
  auto I4d   = Cartesian3d::I4d();
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps / 3.0 * I;
  xt::noalias(Sig) = m_K * treps * I + 2.0 * m_G * Epsd;
  xt::noalias(C) = m_K * II + 2.0 * m_G * I4d;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T2,T4> Elastic::Tangent(const T2& Eps) const
{
  T2 Sig;
  T4 C;
  this->tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
