/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatLinearElastic

================================================================================================= */

#ifndef GMATLINEARELASTIC_CARTESIAN3D_H
#define GMATLINEARELASTIC_CARTESIAN3D_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

// =================================================================================================

namespace GMatLinearElastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

// aliases

using T2 = xt::xtensor_fixed<double, xt::xshape<3,3>>;
using T4 = xt::xtensor_fixed<double, xt::xshape<3,3,3,3>>;

// Tensor operations

template <class T> double trace (const T& A);
template <class T> double ddot22(const T& A, const T& B);

// Unit tensors

T2 I();
T4 II();
T4 I4();
T4 I4rt();
T4 I4s();
T4 I4d();

// -------------------------------------------------------------------------------------------------

// Hydrostatic stress/strain

double hydrostatic(const T2& A);

// Equivalent deviatoric stress/stress

double sigeq(const T2& Sig);
double epseq(const T2& Eps);

// Deviator

void deviator(const T2& A, T2& Ad);
T2   Deviator(const T2& A);

// Matrix version of the functions above (no allocation)

void hydrostatic(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Am);

void sigeq(const xt::xtensor<double,4>& Sig, xt::xtensor<double,2>& sigd);
void epseq(const xt::xtensor<double,4>& Eps, xt::xtensor<double,2>& epsd);

void deviator(const xt::xtensor<double,4>& A, xt::xtensor<double,4>& Ad);

// Auto-allocation allocation of the functions above

xt::xtensor<double,2> Hydrostatic(const xt::xtensor<double,4>& A);

xt::xtensor<double,2> Sigeq(const xt::xtensor<double,4>& Sig);
xt::xtensor<double,2> Epseq(const xt::xtensor<double,4>& Eps);

xt::xtensor<double,4> Deviator(const xt::xtensor<double,4>& Sig);

// -------------------------------------------------------------------------------------------------

class Elastic
{
public:

  // Constructors
  Elastic() = default;
  Elastic(double K, double G);

  // Parameters
  double K() const;
  double G() const;

  // Stress
  void stress(const T2& Eps, T2& Sig) const;

  // Stress & tangent
  void tangent(const T2& Eps, T2& Sig, T4& C) const;

  // Auto-allocation of the functions above
  T2 Stress(const T2& Eps) const;
  std::tuple<T2,T4> Tangent(const T2& Eps) const;

private:

  double m_K; // bulk modulus
  double m_G; // shear modulus
};

// -------------------------------------------------------------------------------------------------

class Matrix
{
public:

  // Constructors

  Matrix() = default;
  Matrix(size_t nelem, size_t nip);
  Matrix(size_t nelem, size_t nip, double K, double G);

  // Shape

  size_t nelem() const;
  size_t nip() const;

  // Parameters

  xt::xtensor<double,2> K() const;
  xt::xtensor<double,2> G() const;

  // Matrix of unit tensors

  xt::xtensor<double,4> I() const;
  xt::xtensor<double,6> II() const;
  xt::xtensor<double,6> I4() const;
  xt::xtensor<double,6> I4rt() const;
  xt::xtensor<double,6> I4s() const;
  xt::xtensor<double,6> I4d() const;

  // Check that a type has been set everywhere (throws if unset points are found)

  void check() const;

  // Set parameters for a batch of points

  void set(const xt::xtensor<size_t,2>& I, double K, double G);

  // Compute (no allocation)

  void stress(
    const xt::xtensor<double,4>& Eps,
          xt::xtensor<double,4>& Sig) const;

  void tangent(
    const xt::xtensor<double,4>& Eps,
          xt::xtensor<double,4>& Sig,
          xt::xtensor<double,6>& C) const;

  // Auto-allocation of the functions above

  xt::xtensor<double,4> Stress(
    const xt::xtensor<double,4>& Eps) const;

  std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Tangent(
    const xt::xtensor<double,4>& Eps) const;

private:

  // Parameters
  xt::xtensor<int   ,2> m_set;
  xt::xtensor<double,2> m_K;
  xt::xtensor<double,2> m_G;

  // Shape
  size_t m_nelem;
  size_t m_nip;
  static const size_t m_ndim=3;
};

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.hpp"

// -------------------------------------------------------------------------------------------------

#endif
