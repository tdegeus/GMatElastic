/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTIC_CARTESIAN3D_H
#define GMATELASTIC_CARTESIAN3D_H

#include <GMatTensor/Cartesian3d.h>

#include "config.h"
#include "version.h"

namespace GMatElastic {

/**
Implementation in a 3-d Cartesian coordinate frame.
*/
namespace Cartesian3d {

/**
Von Mises equivalent strain: norm of strain deviator

\f$ \sqrt{\frac{2}{3} (dev(A))_{ij} (dev(A))_{ji}} \f$

To write to allocated data use epseq().

\param A [..., 3, 3] array.
\return [...] array.
*/
template <class T>
inline auto Epseq(const T& A) -> typename GMatTensor::allocate<xt::get_rank<T>::value - 2, T>::type
{
    return xt::eval(std::sqrt(2.0 / 3.0) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

/**
Same as epseq(), but writes to externally allocated output.

\param A [..., 3, 3] array.
\param ret output [...] array
*/
template <class T, class U>
inline void epseq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(2.0 / 3.0);
}

/**
Von Mises equivalent stress: norm of strain deviator

\f$ \sqrt{\frac{3}{2} (dev(A))_{ij} (dev(A))_{ji}} \f$

To write to allocated data use sigeq().

\param A [..., 3, 3] array.
\return [...] array.
*/
template <class T>
inline auto Sigeq(const T& A) -> typename GMatTensor::allocate<xt::get_rank<T>::value - 2, T>::type
{
    return xt::eval(std::sqrt(1.5) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

/**
Same as Sigeq(), but writes to externally allocated output.

\param A [..., 3, 3] array.
\param ret output [...] array
*/
template <class T, class U>
inline void sigeq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(1.5);
}

/**
Array of material points with a linear elastic constitutive response.
\tparam N Rank of the array.
*/
template <size_t N>
class Elastic : public GMatTensor::Cartesian3d::Array<N> {
protected:
    array_type::tensor<double, N> m_K; ///< Bulk modulus per item.
    array_type::tensor<double, N> m_G; ///< Shear modulus per item.
    array_type::tensor<double, N + 2> m_Eps; ///< Strain tensor per item.
    array_type::tensor<double, N + 2> m_Sig; ///< Stress tensor per item.
    array_type::tensor<double, N + 4> m_C; ///< Tangent per item.

    using GMatTensor::Cartesian3d::Array<N>::m_ndim;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor4;
    using GMatTensor::Cartesian3d::Array<N>::m_size;
    using GMatTensor::Cartesian3d::Array<N>::m_shape;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor4;

public:
    using GMatTensor::Cartesian3d::Array<N>::rank;

    Elastic() = default;

    /**
    Construct system.
    \param K Bulk modulus per item.
    \param G Shear modulus per item.
    */
    template <class T>
    Elastic(const T& K, const T& G)
    {
        GMATELASTIC_ASSERT(K.dimension() == N);
        GMATELASTIC_ASSERT(xt::has_shape(K, G.shape()));
        std::copy(K.shape().cbegin(), K.shape().cend(), m_shape.begin());
        this->init(m_shape);

        m_K = K;
        m_G = G;
        m_Eps = xt::zeros<double>(m_shape_tensor2);
        m_Sig = xt::zeros<double>(m_shape_tensor2);
        m_C = xt::empty<double>(m_shape_tensor4);

#pragma omp parallel
        {
            auto C = xt::adapt(m_C.data(), {m_ndim, m_ndim, m_ndim, m_ndim});
            double K;
            double G;
            auto II = GMatTensor::Cartesian3d::II();
            auto I4d = GMatTensor::Cartesian3d::I4d();

#pragma omp for
            for (size_t i = 0; i < m_size; ++i) {
                C.reset_buffer(&m_C.flat(i * m_stride_tensor4), m_stride_tensor4);
                K = m_K.flat(i);
                G = m_G.flat(i);
                C = K * II + 2.0 * G * I4d;
            }
        }
    }

    /**
    Bulk modulus per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& K() const
    {
        return m_K;
    }

    /**
    Shear modulus per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& G() const
    {
        return m_G;
    }

    /**
    Set strain tensors.
    Internally, this calls refresh() to update stress.
    \tparam T e.g. `array_type::tensor<double, N + 2>`
    \param arg Strain tensor per item [shape(), 3, 3].
    */
    template <class T>
    void set_Eps(const T& arg)
    {
        GMATELASTIC_ASSERT(xt::has_shape(arg, m_shape_tensor2));
        std::copy(arg.cbegin(), arg.cend(), m_Eps.begin());
        this->refresh();
    }

    /**
    Recompute stress from strain.

    From C++, this function need **never** be called: the API takes care of this.

    For Python, this function should **only** be called when you modify elements of Eps().
    For example

        mat.Eps[e, q, 0, 1] = value
        ...
        mat.refresh() # "Eps" was changed without "mat" knowing

    Instead, if you write an nd-array, the API takes care of the refresh. I.e.

        mat.Eps = new_Eps
        # no further action needed, "mat" was refreshed

    Note though that you can call this function as often as you like, you will only loose time.
    */
    virtual void refresh()
    {
#pragma omp parallel for
        for (size_t i = 0; i < m_size; ++i) {

            double K = m_K.flat(i);
            double G = m_G.flat(i);

            const double* Eps = &m_Eps.flat(i * m_stride_tensor2);
            double* Sig = &m_Sig.flat(i * m_stride_tensor2);

            double epsm = GMatTensor::Cartesian3d::pointer::Hydrostatic(Eps);

            Sig[0] = (3.0 * K - 2.0 * G) * epsm + 2.0 * G * Eps[0];
            Sig[1] = 2.0 * G * Eps[1];
            Sig[2] = 2.0 * G * Eps[2];
            Sig[3] = 2.0 * G * Eps[3];
            Sig[4] = (3.0 * K - 2.0 * G) * epsm + 2.0 * G * Eps[4];
            Sig[5] = 2.0 * G * Eps[5];
            Sig[6] = 2.0 * G * Eps[6];
            Sig[7] = 2.0 * G * Eps[7];
            Sig[8] = (3.0 * K - 2.0 * G) * epsm + 2.0 * G * Eps[8];
        }
    }

    /**
    Strain tensor per item.
    \return [shape(), 3, 3].
    */
    const array_type::tensor<double, N + 2>& Eps() const
    {
        return m_Eps;
    }

    /**
    Strain tensor per item.
    The user is responsible for calling refresh() after modifying entries.
    \return [shape(), 3, 3].
    */
    array_type::tensor<double, N + 2>& Eps()
    {
        return m_Eps;
    }

    /**
    Stress tensor per item.
    \return [shape(), 3, 3].
    */
    const array_type::tensor<double, N + 2>& Sig() const
    {
        return m_Sig;
    }

    /**
    Tangent tensor per item.
    \return [shape(), 3, 3, 3, 3].
    */
    const array_type::tensor<double, N + 4>& C() const
    {
        return m_C;
    }

    /**
    Potential energy per item.
    \return [shape()].
    */
    virtual array_type::tensor<double, N> energy() const
    {
        array_type::tensor<double, N> ret = xt::empty<double>(m_shape);
        namespace GT = GMatTensor::Cartesian3d::pointer;

#pragma omp parallel for
        for (size_t i = 0; i < m_size; ++i) {

            double K = m_K.flat(i);
            double G = m_G.flat(i);

            const double* Eps = &m_Eps.flat(i * m_stride_tensor2);

            std::array<double, m_stride_tensor2> Epsd;
            double epsm = GT::Hydrostatic_deviatoric(Eps, &Epsd[0]);
            double epsd = std::sqrt(0.5 * GT::A2s_ddot_B2s(&Epsd[0], &Epsd[0]));

            double U = 3.0 * K * std::pow(epsm, 2.0);
            double V = 2.0 * G * std::pow(epsd, 2.0);

            ret.flat(i) = U + V;
        }

        return ret;
    }
};

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
