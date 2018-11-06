/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatLinearElastic

================================================================================================= */

#ifndef GMATLINEARELASTIC_CARTESIAN3D_H
#define GMATLINEARELASTIC_CARTESIAN3D_H

// -------------------------------------------------------------------------------------------------

#include "GMatLinearElastic.h"

// =================================================================================================

namespace GMatLinearElastic {

// -------------------------------------------------------------------------------------------------

class Cartesian3d
{
public:

  // constructors
  Cartesian3d() = default;
  Cartesian3d(size_t nelem, size_t nip);
  Cartesian3d(size_t nelem, size_t nip, double kappa, double mu);

  // return shape
  size_t nelem() const;
  size_t nip() const;

  // parameters
  xt::xtensor<double,2> kappa() const;
  xt::xtensor<double,2> mu() const;

  // check that a type has been set everywhere
  void check() const;

  // set parameters
  void set(const xt::xtensor<size_t,2> &I, double kappa, double mu);

  // compute (no allocation)
  void Sig      (const xt::xtensor<double,4> &a_Eps  , xt::xtensor<double,4> &a_Sig      ) const;
  void Sig2d    (const xt::xtensor<double,4> &a_Eps2d, xt::xtensor<double,4> &a_Sig2d    ) const;
  void Tangent  (                                      xt::xtensor<double,6> &a_Tangent  ) const;
  void Tangent2d(                                      xt::xtensor<double,6> &a_Tangent2d) const;

  // compute (return allocated result)
  xt::xtensor<double,4> Sig      (const xt::xtensor<double,4> &a_Eps  ) const;
  xt::xtensor<double,4> Sig2d    (const xt::xtensor<double,4> &a_Eps2d) const;
  xt::xtensor<double,6> Tangent  (                                    ) const;
  xt::xtensor<double,6> Tangent2d(                                    ) const;

private:

  xt::xtensor<int   ,2> m_set;
  xt::xtensor<double,2> m_kappa;
  xt::xtensor<double,2> m_mu;
  size_t                m_nelem;
  size_t                m_nip;
};

// -------------------------------------------------------------------------------------------------

namespace Support {
namespace Cartesian2d {

using T2 = xt::xtensor_fixed<double, xt::xshape<2,2>>;
using T4 = xt::xtensor_fixed<double, xt::xshape<2,2,2,2>>;

template<class T> double trace (const T &A);
template<class T> double ddot22(const T &A, const T &B);

T2 I();
T4 II();
T4 I4();
T4 I4rt();
T4 I4s();
T4 I4d();

}}

// -------------------------------------------------------------------------------------------------

namespace Support {
namespace Cartesian3d {

using T2 = xt::xtensor_fixed<double, xt::xshape<3,3>>;
using T4 = xt::xtensor_fixed<double, xt::xshape<3,3,3,3>>;

template<class T> double trace (const T &A);
template<class T> double ddot22(const T &A, const T &B);

T2 I();
T4 II();
T4 I4();
T4 I4rt();
T4 I4s();
T4 I4d();

}}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
