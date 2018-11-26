/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatLinearElastic

================================================================================================= */

#ifndef GMATLINEARELASTIC_CARTESIAN3D_HPP
#define GMATLINEARELASTIC_CARTESIAN3D_HPP

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.h"

// =================================================================================================

namespace GMatLinearElastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

template<class T>
inline double trace(const T &A)
{
  return A(0,0) + A(1,1) + A(2,2);
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline double ddot22(const T &A, const T &B)
{
  return A(0,0) * B(0,0) + 2.0 * A(0,1) * B(0,1) + 2.0 * A(0,2) * B(0,2) +
         A(1,1) * B(1,1) + 2.0 * A(1,2) * B(1,2) +
         A(2,2) * B(2,2);
}

// -------------------------------------------------------------------------------------------------

inline T2 I()
{
  return T2({{1., 0., 0.},
             {0., 1., 0.},
             {0., 0., 1.}});
}

// -------------------------------------------------------------------------------------------------

inline T4 II()
{
  T4 out;

  out.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          if ( i == j and k == l )
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4()
{
  T4 out;

  out.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          if ( i == l and j == k )
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4rt()
{
  T4 out;

  out.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          if ( i == k and j == l )
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4s()
{
  return .5 * ( I4() + I4rt() );
}

// -------------------------------------------------------------------------------------------------

inline T4 I4d()
{
  return I4s() - II()/3.;
}

// -------------------------------------------------------------------------------------------------

inline Elastic::Elastic(double kappa, double mu) : m_kappa(kappa), m_mu(mu) { }

inline double Elastic::kappa() const { return m_kappa; }

inline double Elastic::mu() const { return m_mu; }

// -------------------------------------------------------------------------------------------------

inline T2 Elastic::Sig(const T2 &Eps) const
{
  // define identity tensor
  T2 I = Cartesian3d::I();

  // decompose strain
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps/3. * I;

  // compute stress
  return m_kappa * treps * I + 2. * m_mu * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T2,T4> Elastic::Tangent(const T2 &Eps) const
{
  // define identity tensor
  T2 I   = Cartesian3d::I();
  T4 II  = Cartesian3d::II();
  T4 I4d = Cartesian3d::I4d();

  // decompose strain
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps/3. * I;

  // compute stress
  auto Sig = m_kappa * treps * I + 2. * m_mu * Epsd;

  // compute tangent
  auto C4 = m_kappa * II + 2. * m_mu * I4d;

  return std::make_tuple(Sig, C4);
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::nelem() const { return m_nelem; }

inline size_t Matrix::nip() const { return m_nip; }

inline xt::xtensor<double,2> Matrix::kappa() const { return m_kappa; };

inline xt::xtensor<double,2> Matrix::mu() const { return m_mu; };

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip) : m_nelem(nelem), m_nip(nip)
{
  m_set   = xt::zeros<int   >({nelem, nip});
  m_kappa = xt::empty<double>({nelem, nip});
  m_mu    = xt::empty<double>({nelem, nip});
}

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip, double kappa, double mu) :
  m_nelem(nelem), m_nip(nip)
{
  m_set   = xt::ones<int   >({nelem, nip});
  m_kappa = xt::ones<double>({nelem, nip}) * kappa;
  m_mu    = xt::ones<double>({nelem, nip}) * mu;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::check() const
{
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t q = 0 ; q < m_nip ; ++q )
      if ( m_set(e,q) == 0 )
        throw std::runtime_error("No type set for: "+std::to_string(e)+", "+std::to_string(q));
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::set(const xt::xtensor<size_t,2> &I, double kappa, double mu)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_set.shape() );
    // - empty material definitions
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t q = 0 ; q < m_nip ; ++q )
        if ( I(e,q) )
          assert( m_set(e,q) == 0 );
  #endif

  // set type and position in material vector
  for ( size_t e = 0 ; e < m_nelem ; ++e ) {
    for ( size_t q = 0 ; q < m_nip ; ++q ) {
      if ( I(e,q) ) {
        m_set  (e,q) = 1;
        m_kappa(e,q) = kappa;
        m_mu   (e,q) = mu;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::I() const
{
  xt::xtensor<double,4> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T2 unit = Cartesian3d::I();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0), xt::xshape<m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::II() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::II();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4rt() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4rt();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4s() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4s();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4d() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4d();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::Sig(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Sig) const
{
  assert( a_Eps.shape()[0] == m_nelem       );
  assert( a_Eps.shape()[1] == m_nip         );
  assert( a_Eps.shape()[2] == m_ndim        );
  assert( a_Eps.shape()[3] == m_ndim        );
  assert( a_Eps.shape()    == a_Sig.shape() );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensor
    T2 I = Cartesian3d::I();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        // alias
        auto  Eps   = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  Sig   = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto& kappa = m_kappa(e,q);
        auto& mu    = m_mu   (e,q);

        // decompose strain
        auto  treps = trace(Eps);
        auto  Epsd  = Eps - treps/3. * I;

        // compute stress
        xt::noalias(Sig) = kappa * treps * I + 2. * mu * Epsd;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::Tangent(const xt::xtensor<double,4> &a_Eps,
  xt::xtensor<double,4> &a_Sig, xt::xtensor<double,6> &a_Tangent) const
{
  assert( a_Eps.shape()[0]     == m_nelem       );
  assert( a_Eps.shape()[1]     == m_nip         );
  assert( a_Eps.shape()[2]     == m_ndim        );
  assert( a_Eps.shape()[3]     == m_ndim        );
  assert( a_Eps.shape()        == a_Sig.shape() );
  assert( a_Tangent.shape()[0] == m_nelem       );
  assert( a_Tangent.shape()[1] == m_nip         );
  assert( a_Tangent.shape()[2] == m_ndim        );
  assert( a_Tangent.shape()[3] == m_ndim        );
  assert( a_Tangent.shape()[4] == m_ndim        );
  assert( a_Tangent.shape()[5] == m_ndim        );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensors
    T2 I   = Cartesian3d::I();
    T4 II  = Cartesian3d::II();
    T4 I4d = Cartesian3d::I4d();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        // alias
        auto  Eps   = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  Sig   = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  C4    = xt::adapt(&a_Tangent(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        auto& kappa = m_kappa(e,q);
        auto& mu    = m_mu   (e,q);

        // decompose strain
        auto  treps = trace(Eps);
        auto  Epsd  = Eps - treps/3. * I;

        // compute stress
        xt::noalias(Sig) = kappa * treps * I + 2. * mu * Epsd;

        // compute tangent
        xt::noalias(C4) = kappa * II + 2. * mu * I4d;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::Sig(const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,4> a_Sig = xt::empty<double>(a_Eps.shape());

  this->Sig(a_Eps, a_Sig);

  return a_Sig;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Matrix::Tangent(
  const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,4> a_Sig     = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim});
  xt::xtensor<double,6> a_Tangent = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim,m_ndim,m_ndim});

  this->Tangent(a_Eps, a_Sig, a_Tangent);

  return std::make_tuple(a_Sig, a_Tangent);
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif
