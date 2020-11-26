/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_H
#define GMATELASTIC_CARTESIAN3D_H

#include <GMatTensor/Cartesian3d.h>

#include "config.h"

namespace GMatElastic {
namespace Cartesian3d {

// Unit tensors

using GMatTensor::Cartesian3d::I2;
using GMatTensor::Cartesian3d::II;
using GMatTensor::Cartesian3d::I4;
using GMatTensor::Cartesian3d::I4rt;
using GMatTensor::Cartesian3d::I4s;
using GMatTensor::Cartesian3d::I4d;

// Tensor decomposition

using GMatTensor::Cartesian3d::hydrostatic;
using GMatTensor::Cartesian3d::Hydrostatic;
using GMatTensor::Cartesian3d::deviatoric;
using GMatTensor::Cartesian3d::Deviatoric;

// Equivalent strain

template <class T, class U>
inline void epseq(const T& A, U& B);

template <class T>
inline auto Epseq(const T& A);

// Equivalent stress

template <class T, class U>
inline void sigeq(const T& A, U& B);

template <class T>
inline auto Sigeq(const T& A);

// Material point

class Elastic
{
public:
    // Constructors
    Elastic() = default;
    Elastic(double K, double G);

    // Parameters
    double K() const;
    double G() const;

    // Set strain
    template <class T>
    void setStrain(const T& Eps);

    template <class T>
    void setStrainIterator(const T& begin); // presumes: contiguous + row-major & symmetric

    // Stress (no allocation, overwrites "Sig" / writes to "begin")
    template <class T>
    void stress(T& Sig) const;

    template <class T>
    void stressIterator(const T& begin) const; // presumes: contiguous + row-major

    // Tangent (no allocation, overwrites "C")
    template <class T>
    void tangent(T& C) const;

    // Auto-allocation
    xt::xtensor<double, 2> Stress() const;
    xt::xtensor<double, 4> Tangent() const;

private:
    double m_K;                  // bulk modulus
    double m_G;                  // shear modulus
    std::array<double, 9> m_Eps; // strain tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    std::array<double, 9> m_Sig; // stress tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
};

// Material identifier

struct Type {
    enum Value {
        Unset,
        Elastic,
    };
};

// Array of material points

template <size_t N>
class Array : public GMatTensor::Cartesian3d::Array<N>
{
public:
    using GMatTensor::Cartesian3d::Array<N>::rank;

    // Constructors

    Array() = default;
    Array(const std::array<size_t, N>& shape);
    Array(const std::array<size_t, N>& shape, double K, double G);

    // Overloaded methods

    /*
    std::array<size_t, N> shape() const;

    xt::xtensor<double, N + 2> I2() const;
    xt::xtensor<double, N + 4> II() const;
    xt::xtensor<double, N + 4> I4() const;
    xt::xtensor<double, N + 4> I4rt() const;
    xt::xtensor<double, N + 4> I4s() const;
    xt::xtensor<double, N + 4> I4d() const;
    */

    // Type

    xt::xtensor<size_t, N> type() const;
    xt::xtensor<size_t, N> isElastic() const;

    // Parameters

    xt::xtensor<double, N> K() const;
    xt::xtensor<double, N> G() const;

    // Check that a type has been set everywhere (throws if unset points are found)

    void check() const;

    // Set parameters for a batch of points

    void setElastic(const xt::xtensor<size_t, N>& I, double K, double G);

    // Set strain tensor, get the response

    void setStrain(const xt::xtensor<double, N + 2>& Eps);
    void stress(xt::xtensor<double, N + 2>& Sig) const;
    void tangent(xt::xtensor<double, N + 4>& C) const;

    // Auto-allocation of the functions above

    xt::xtensor<double, N + 2> Stress() const;
    xt::xtensor<double, N + 4> Tangent() const;

private:
    // Material vectors
    std::vector<Elastic> m_Elastic;

    // Identifiers for each matrix entry
    xt::xtensor<size_t, N> m_type;  // type (e.g. "Type::Elastic")
    xt::xtensor<size_t, N> m_index; // index from the relevant material vector (e.g. "m_Elastic")

    // Shape
    using GMatTensor::Cartesian3d::Array<N>::m_ndim;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor4;
    using GMatTensor::Cartesian3d::Array<N>::m_size;
    using GMatTensor::Cartesian3d::Array<N>::m_shape;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor4;

    // Internal check
    bool m_allSet = false; // true if all points have a material assigned
    void checkAllSet(); // check if all points have a material assigned (modifies "m_allSet")
};

} // namespace Cartesian3d
} // namespace GMatElastic

#include "Cartesian3d.hpp"
#include "Cartesian3d_Array.hpp"
#include "Cartesian3d_Elastic.hpp"

#endif
