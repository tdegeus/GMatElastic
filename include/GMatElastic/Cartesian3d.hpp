/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_HPP
#define GMATELASTIC_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastic {
namespace Cartesian3d {

namespace detail {
namespace xtensor {

    template <class T>
    inline auto trace(const T& A)
    {
        return A(0, 0) + A(1, 1) + A(2, 2);
    }

    template <class T, class U>
    inline auto A2_ddot_B2(const T& A, const U& B)
    {
        return A(0, 0) * B(0, 0)
             + A(1, 1) * B(1, 1)
             + A(2, 2) * B(2, 2)
             + 2.0 * A(0, 1) * B(0, 1)
             + 2.0 * A(0, 2) * B(0, 2)
             + 2.0 * A(1, 2) * B(1, 2);
    }

} // namespace xtensor
} // namespace detail

namespace detail {
namespace pointer {

    template <class T>
    inline auto trace(const T A)
    {
        return A[0] + A[4] + A[8];
    }

    template <class T, class U>
    inline void deviatoric(const T A, U ret)
    {
        auto m = 0.5 * (A[0] + A[4] + A[8]);
        ret[0] = A[0] - m;
        ret[1] = A[1];
        ret[2] = A[2];
        ret[3] = A[3];
        ret[4] = A[4] - m;
        ret[5] = A[5];
        ret[6] = A[6];
        ret[7] = A[7];
        ret[8] = A[8] - m;
    }

    template <class T>
    inline auto deviatoric_ddot_deviatoric(const T A)
    {
        auto m = 0.5 * (A[0] + A[4] + A[8]);
        return (A[0] - m) * (A[0] - m)
             + (A[4] - m) * (A[4] - m)
             + (A[8] - m) * (A[8] - m)
             + 2.0 * A[1] * A[1]
             + 2.0 * A[2] * A[2]
             + 2.0 * A[5] * A[5];
    }

} // namespace pointer
} // namespace detail

namespace detail {

    template <class T>
    inline T trace(const std::array<T, 9>& A)
    {
        return A[0] + A[4] + A[8];
    }

    template <class T>
    inline T hydrostatic_deviator(const std::array<T, 9>& A, std::array<T, 9>& ret)
    {
        auto m = 0.5 * (A[0] + A[4] + A[8]);
        ret[0] = A[0] - m;
        ret[1] = A[1];
        ret[2] = A[2];
        ret[3] = A[3];
        ret[4] = A[4] - m;
        ret[5] = A[5];
        ret[6] = A[6];
        ret[7] = A[7];
        ret[8] = A[8] - m;
        return m;
    }

    template <class T>
    inline T A2_ddot_B2(const std::array<T, 9>& A, const std::array<T, 9>& B)
    {
        auto m = 0.5 * (A[0] + A[4] + A[8]);
        return (A[0] - m) * (B[0] - m)
             + (A[4] - m) * (B[4] - m)
             + (A[8] - m) * (B[8] - m)
             + 2.0 * A[1] * B[1]
             + 2.0 * A[2] * B[2]
             + 2.0 * A[5] * B[5];
    }

} // namespace detail

inline Tensor2 I2()
{
    return Tensor2({{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}});
}

inline Tensor4 II()
{
    Tensor4 ret = Tensor4::from_shape({3, 3, 3, 3});
    ret.fill(0.0);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    if (i == j && k == l) {
                        ret(i, j, k, l) = 1.0;
                    }
                }
            }
        }
    }

    return ret;
}

inline Tensor4 I4()
{
    Tensor4 ret = Tensor4::from_shape({3, 3, 3, 3});
    ret.fill(0.0);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    if (i == l && j == k) {
                        ret(i, j, k, l) = 1.0;
                    }
                }
            }
        }
    }

    return ret;
}

inline Tensor4 I4rt()
{
    Tensor4 ret = Tensor4::from_shape({3, 3, 3, 3});
    ret.fill(0.0);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    if (i == k && j == l) {
                        ret(i, j, k, l) = 1.0;
                    }
                }
            }
        }
    }

    return ret;
}

inline Tensor4 I4s()
{
    return 0.5 * (I4() + I4rt());
}

inline Tensor4 I4d()
{
    return I4s() - II() / 3.0;
}

namespace detail {

    template <class T>
    struct equiv_impl
    {
        using value_type = typename T::value_type;
        using shape_type = typename T::shape_type;
        static_assert(xt::has_fixed_rank_t<T>::value, "Only fixed rank allowed.");
        static_assert(xt::get_rank<T>::value >= 2, "Rank too low.");
        static const size_t rank = xt::get_rank<T>::value;
        static const size_t scalar_rank = rank - 2;
        static const size_t ndim = 3;
        static const size_t stride = ndim * ndim;

        template <class S>
        static size_t getMatrixSize(const S& arg)
        {
            size_t ret = 1;
            for (size_t i = 0; i < scalar_rank; ++i) {
                ret *= arg[i];
            }
            return ret;
        }

        template <class S>
        static std::array<size_t, rank> getShape(const S& arg)
        {
            std::array<size_t, rank> ret;
            for (size_t i = 0; i < rank; ++i) {
                ret[i] = arg[i];
            }
            return ret;
        }

        template <class S>
        static std::array<size_t, scalar_rank> getShapeScalar(const S& arg)
        {
            std::array<size_t, scalar_rank> ret;
            for (size_t i = 0; i < scalar_rank; ++i) {
                ret[i] = arg[i];
            }
            return ret;
        }

        template <class S>
        static std::array<size_t, rank> getShapeTensor(const S& arg)
        {
            std::array<size_t, rank> ret;
            for (size_t i = 0; i < scalar_rank; ++i) {
                ret[i] = arg[i];
            }
            for (size_t i = scalar_rank; i < rank; ++i) {
                ret[i] = ndim;
            }
            return ret;
        }

        static void deviatoric_no_alloc(const T& A, xt::xtensor<value_type, rank>& B)
        {
            GMATELASTOPLASTICQPOT_ASSERT(getShape(A.shape()) == getShapeTensor(B.shape()));
            GMATELASTOPLASTICQPOT_ASSERT(getShape(A.shape()) == getShape(B.shape()));
            #pragma omp parallel for
            for (size_t i = 0; i < getMatrixSize(A.shape()); ++i) {
                detail::pointer::deviatoric(&A.data()[i * stride], &B.data()[i * stride]);
            }
        }

        static void hydrostatic_no_alloc(const T& A, xt::xtensor<value_type, scalar_rank>& B)
        {
            GMATELASTOPLASTICQPOT_ASSERT(getShape(A.shape()) == getShapeTensor(B.shape()));
            #pragma omp parallel for
            for (size_t i = 0; i < getMatrixSize(A.shape()); ++i) {
                B.data()[i] = detail::pointer::trace(&A.data()[i * stride]) / 3.0;
            }
        }

        static void epseq_no_alloc(const T& A, xt::xtensor<value_type, scalar_rank>& B)
        {
            GMATELASTOPLASTICQPOT_ASSERT(getShape(A.shape()) == getShapeTensor(B.shape()));
            #pragma omp parallel for
            for (size_t i = 0; i < getMatrixSize(A.shape()); ++i) {
                auto b = detail::pointer::deviatoric_ddot_deviatoric(&A.data()[i * stride]);
                B.data()[i] = std::sqrt(2.0 / 3.0 * b);
            }
        }

        static void sigeq_no_alloc(const T& A, xt::xtensor<value_type, scalar_rank>& B)
        {
            GMATELASTOPLASTICQPOT_ASSERT(getShape(A.shape()) == getShapeTensor(B.shape()));
            #pragma omp parallel for
            for (size_t i = 0; i < getMatrixSize(A.shape()); ++i) {
                auto b = detail::pointer::deviatoric_ddot_deviatoric(&A.data()[i * stride]);
                B.data()[i] = std::sqrt(1.5 * b);
            }
        }

        static auto deviatoric_alloc(const T& A)
        {
            xt::xtensor<value_type, rank> B = xt::empty<value_type>(A.shape());
            deviatoric_no_alloc(A, B);
            return B;
        }

        static auto hydrostatic_alloc(const T& A)
        {
            xt::xtensor<value_type, scalar_rank> B = xt::empty<value_type>(getShapeScalar(A.shape()));
            hydrostatic_no_alloc(A, B);
            return B;
        }

        static auto epseq_alloc(const T& A)
        {
            xt::xtensor<value_type, scalar_rank> B = xt::empty<value_type>(getShapeScalar(A.shape()));
            epseq_no_alloc(A, B);
            return B;
        }

        static auto sigeq_alloc(const T& A)
        {
            xt::xtensor<value_type, scalar_rank> B = xt::empty<value_type>(getShapeScalar(A.shape()));
            sigeq_no_alloc(A, B);
            return B;
        }
    };

} // namespace detail

template <class T, class U>
inline void hydrostatic(const T& A, U& B)
{
    return detail::equiv_impl<T>::hydrostatic_no_alloc(A, B);
}

template <class T>
inline auto Hydrostatic(const T& A)
{
    return detail::equiv_impl<T>::hydrostatic_alloc(A);
}

template <class T, class U>
inline void deviatoric(const T& A, U& B)
{
    return detail::equiv_impl<T>::deviatoric_no_alloc(A, B);
}

template <class T>
inline auto Deviatoric(const T& A)
{
    return detail::equiv_impl<T>::deviatoric_alloc(A);
}

template <class T, class U>
inline void epseq(const T& A, U& B)
{
    return detail::equiv_impl<T>::epseq_no_alloc(A, B);
}

template <class T>
inline auto Epseq(const T& A)
{
    return detail::equiv_impl<T>::epseq_alloc(A);
}

template <class T, class U>
inline void sigeq(const T& A, U& B)
{
    return detail::equiv_impl<T>::sigeq_no_alloc(A, B);
}

template <class T>
inline auto Sigeq(const T& A)
{
    return detail::equiv_impl<T>::sigeq_alloc(A);
}

} // namespace Cartesian3d
} // namespace GMatElastic

#endif
