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
inline double trace(const T& A)
{
  return A(0,0) + A(1,1) + A(2,2);
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline double ddot22(const T& A, const T& B)
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

  for (size_t i = 0 ; i < 3 ; ++i)
    for (size_t j = 0 ; j < 3 ; ++j)
      for (size_t k = 0 ; k < 3 ; ++k)
        for (size_t l = 0 ; l < 3 ; ++l)
          if (i == j and k == l)
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4()
{
  T4 out;

  out.fill(0.0);

  for (size_t i = 0 ; i < 3 ; ++i)
    for (size_t j = 0 ; j < 3 ; ++j)
      for (size_t k = 0 ; k < 3 ; ++k)
        for (size_t l = 0 ; l < 3 ; ++l)
          if (i == l and j == k)
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4rt()
{
  T4 out;

  out.fill(0.0);

  for (size_t i = 0 ; i < 3 ; ++i)
    for (size_t j = 0 ; j < 3 ; ++j)
      for (size_t k = 0 ; k < 3 ; ++k)
        for (size_t l = 0 ; l < 3 ; ++l)
          if (i == k and j == l)
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

inline double epseq(const T2& Eps)
{
  T2 Epsd = Eps - trace(Eps)/3. * I();
  return std::sqrt(2./3.*ddot22(Epsd,Epsd));
}

// -------------------------------------------------------------------------------------------------

inline double sigeq(const T2& Sig)
{
  T2 Sigd = Sig - trace(Sig)/3. * I();
  return std::sqrt(1.5*ddot22(Sigd,Sigd));
}

// -------------------------------------------------------------------------------------------------

inline void deviator(const T2& A, T2 &Ad)
{
  xt::noalias(Ad) = A - trace(A)/3. * I();
}

// -------------------------------------------------------------------------------------------------

inline T2 Deviator(const T2& A)
{
  return A - trace(A)/3. * I();
}

// -------------------------------------------------------------------------------------------------

inline double hydrostatic(const T2& A)
{
  return trace(A)/3.;
}

// -------------------------------------------------------------------------------------------------

inline void epseq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq)
{
  assert(A.shape()[0] == Aeq.shape()[0]);
  assert(A.shape()[1] == Aeq.shape()[1]);
  assert(A.shape()[2] == 3);
  assert(A.shape()[3] == 3);

  #pragma omp parallel
  {
    T2 I = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0 ; e < A.shape()[0] ; ++e) {
      for (size_t q = 0 ; q < A.shape()[1] ; ++q) {
        auto Ai  = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        auto Aid = Ai - trace(Ai)/3. * I;
        Aeq(e,q) = std::sqrt(2./3.*ddot22(Aid,Aid));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void sigeq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq)
{
  assert(A.shape()[0] == Aeq.shape()[0]);
  assert(A.shape()[1] == Aeq.shape()[1]);
  assert(A.shape()[2] == 3);
  assert(A.shape()[3] == 3);

  #pragma omp parallel
  {
    T2 I = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0 ; e < A.shape()[0] ; ++e) {
      for (size_t q = 0 ; q < A.shape()[1] ; ++q) {
        auto Ai  = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        auto Aid = Ai - trace(Ai)/3. * I;
        Aeq(e,q) = std::sqrt(1.5*ddot22(Aid,Aid));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void deviator(const xt::xtensor<double,4>& A, xt::xtensor<double,4>& Ad)
{
  assert(A.shape() == Ad.shape());
  assert(A.shape()[2] == 3);
  assert(A.shape()[3] == 3);

  #pragma omp parallel
  {
    T2 I = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0 ; e < A.shape()[0] ; ++e) {
      for (size_t q = 0 ; q < A.shape()[1] ; ++q) {
        auto Ai  = xt::adapt(&A (e,q,0,0), xt::xshape<3,3>());
        auto Aid = xt::adapt(&Ad(e,q,0,0), xt::xshape<3,3>());
        xt::noalias(Aid) = Ai - trace(Ai)/3. * I;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void hydrostatic(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Am)
{
  assert(A.shape()[0] == Am.shape()[0]);
  assert(A.shape()[1] == Am.shape()[1]);
  assert(A.shape()[2] == 3);
  assert(A.shape()[3] == 3);

  #pragma omp parallel
  {
    #pragma omp for
    for (size_t e = 0 ; e < A.shape()[0] ; ++e) {
      for (size_t q = 0 ; q < A.shape()[1] ; ++q) {
        auto Ai = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        Am(e,q) = trace(Ai)/3.;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Epseq(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Aeq = xt::empty<double>({A.shape()[0], A.shape()[1]});
  epseq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Sigeq(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Aeq = xt::empty<double>({A.shape()[0], A.shape()[1]});
  sigeq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Deviator(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,4> Ad = xt::empty<double>(A.shape());
  deviator(A, Ad);
  return Ad;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Hydrostatic(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Am = xt::empty<double>({A.shape()[0], A.shape()[1]});
  hydrostatic(A, Am);
  return Am;
}

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

inline void Elastic::stress(const T2& Eps, T2& Sig) const
{
  T2 I = Cartesian3d::I();
  auto treps = trace(Eps);
  auto Epsd = Eps - treps/3. * I;
  xt::noalias(Sig) = m_K * treps * I + 2. * m_G * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline void Elastic::tangent(const T2& Eps, T2& Sig, T4& C) const
{
  T2 I   = Cartesian3d::I();
  T4 II  = Cartesian3d::II();
  T4 I4d = Cartesian3d::I4d();
  auto treps = trace(Eps);
  auto Epsd = Eps - treps/3. * I;
  xt::noalias(Sig) = m_K * treps * I + 2. * m_G * Epsd;
  xt::noalias(C) = m_K * II + 2. * m_G * I4d;
}

// -------------------------------------------------------------------------------------------------

inline T2 Elastic::Stress(const T2& Eps) const
{
  T2 Sig;
  stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T2,T4> Elastic::Tangent(const T2& Eps) const
{
  T2 Sig;
  T4 C;
  tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip) : m_nelem(nelem), m_nip(nip)
{
  m_set = xt::zeros<int   >({nelem, nip});
  m_K   = xt::empty<double>({nelem, nip});
  m_G   = xt::empty<double>({nelem, nip});
}

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip, double K, double G) : m_nelem(nelem), m_nip(nip)
{
  m_set = xt::ones<int   >({nelem, nip});
  m_K   = xt::ones<double>({nelem, nip}) * K;
  m_G   = xt::ones<double>({nelem, nip}) * G;
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::nip() const
{
  return m_nip;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Matrix::K() const
{
  return m_K;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Matrix::G() const
{
  return m_G;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::I() const
{
  xt::xtensor<double,4> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T2 unit = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
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
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
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
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
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
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
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
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
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
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::check() const
{
  for (size_t e = 0 ; e < m_nelem ; ++e)
    for (size_t q = 0 ; q < m_nip ; ++q)
      if (m_set(e,q) == 0)
        throw std::runtime_error("No type set for: "+std::to_string(e)+", "+std::to_string(q));
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::set(const xt::xtensor<size_t,2>& I, double K, double G)
{
  #ifndef NDEBUG
    assert(I.shape() == m_set.shape());
    for (size_t e = 0 ; e < m_nelem ; ++e)
      for (size_t q = 0 ; q < m_nip ; ++q)
        if (I(e,q))
          assert(m_set(e,q) == 0);
  #endif

  for (size_t e = 0 ; e < m_nelem ; ++e) {
    for (size_t q = 0 ; q < m_nip ; ++q) {
      if (I(e,q)) {
        m_set(e,q) = 1;
        m_K  (e,q) = K;
        m_G  (e,q) = G;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::stress(const xt::xtensor<double,4>& a_Eps, xt::xtensor<double,4>& a_Sig) const
{
  assert(a_Eps.shape()[0] == m_nelem);
  assert(a_Eps.shape()[1] == m_nip);
  assert(a_Eps.shape()[2] == m_ndim);
  assert(a_Eps.shape()[3] == m_ndim);
  assert(a_Eps.shape() == a_Sig.shape());

  #pragma omp parallel
  {
    T2 I = Cartesian3d::I();

    #pragma omp for
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
        auto  Eps = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  Sig = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto& K = m_K(e,q);
        auto& G = m_G(e,q);
        auto  treps = trace(Eps);
        auto  Epsd = Eps - treps/3. * I;
        xt::noalias(Sig) = K * treps * I + 2. * G * Epsd;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::tangent(
  const xt::xtensor<double,4>& a_Eps,
        xt::xtensor<double,4>& a_Sig,
        xt::xtensor<double,6>& a_Tangent) const
{
  assert(a_Eps.shape()[0] == m_nelem);
  assert(a_Eps.shape()[1] == m_nip);
  assert(a_Eps.shape()[2] == m_ndim);
  assert(a_Eps.shape()[3] == m_ndim);
  assert(a_Eps.shape() == a_Sig.shape());
  assert(a_Tangent.shape()[0] == m_nelem);
  assert(a_Tangent.shape()[1] == m_nip);
  assert(a_Tangent.shape()[2] == m_ndim);
  assert(a_Tangent.shape()[3] == m_ndim);
  assert(a_Tangent.shape()[4] == m_ndim);
  assert(a_Tangent.shape()[5] == m_ndim);

  #pragma omp parallel
  {
    T2 I   = Cartesian3d::I();
    T4 II  = Cartesian3d::II();
    T4 I4d = Cartesian3d::I4d();

    #pragma omp for
    for (size_t e = 0 ; e < m_nelem ; ++e) {
      for (size_t q = 0 ; q < m_nip ; ++q) {
        auto  Eps = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  Sig = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  C4 = xt::adapt(&a_Tangent(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        auto& K = m_K(e,q);
        auto& G = m_G(e,q);
        auto  treps = trace(Eps);
        auto  Epsd = Eps - treps/3. * I;
        xt::noalias(Sig) = K * treps * I + 2. * G * Epsd;
        xt::noalias(C4) = K * II + 2. * G * I4d;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::Stress(const xt::xtensor<double,4>& a_Eps) const
{
  xt::xtensor<double,4> a_Sig = xt::empty<double>(a_Eps.shape());
  stress(a_Eps, a_Sig);
  return a_Sig;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Matrix::Tangent(
  const xt::xtensor<double,4>& a_Eps) const
{
  xt::xtensor<double,4> a_Sig     = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim});
  xt::xtensor<double,6> a_Tangent = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim,m_ndim,m_ndim});
  tangent(a_Eps, a_Sig, a_Tangent);
  return std::make_tuple(a_Sig, a_Tangent);
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif
