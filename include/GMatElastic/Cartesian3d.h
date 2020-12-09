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

using GMatTensor::Cartesian3d::O2;
using GMatTensor::Cartesian3d::O4;
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
inline void epseq(const T& A, U& ret);

template <class T>
inline auto Epseq(const T& A);

// Equivalent stress

template <class T, class U>
inline void sigeq(const T& A, U& ret);

template <class T>
inline auto Sigeq(const T& A);

// Material point

class Elastic
{
public:
    Elastic() = default;
    Elastic(double K, double G);

    double K() const;
    double G() const;
    double energy() const;

    template <class T> void setStrain(const T& arg);
    template <class T> void strain(T& ret) const;
    template <class T> void stress(T& ret) const;
    template <class T> void tangent(T& ret) const;

    template <class T> void setStrainPtr(const T* arg);
    template <class T> void strainPtr(T* ret) const;
    template <class T> void stressPtr(T* ret) const;
    template <class T> void tangentPtr(T* ret) const;

    xt::xtensor<double, 2> Strain() const;
    xt::xtensor<double, 2> Stress() const;
    xt::xtensor<double, 4> Tangent() const;

private:
    double m_K; // bulk modulus
    double m_G; // shear modulus
    std::array<double, 9> m_Eps; // strain tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    std::array<double, 9> m_Sig; // stress tensor ,,
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

    // Overloaded methods:
    // - "shape"
    // - unit tensors: "I2", "II", "I4", "I4rt", "I4s", "I4d"

    // Type

    xt::xtensor<size_t, N> type() const;
    xt::xtensor<size_t, N> isElastic() const;

    // Parameters

    xt::xtensor<double, N> K() const;
    xt::xtensor<double, N> G() const;

    // Set parameters for a batch of points

    void setElastic(const xt::xtensor<size_t, N>& I, double K, double G);

    // Set strain tensor, get the response

    void setStrain(const xt::xtensor<double, N + 2>& arg);
    void strain(xt::xtensor<double, N + 2>& ret) const;
    void stress(xt::xtensor<double, N + 2>& ret) const;
    void tangent(xt::xtensor<double, N + 4>& ret) const;

    // Auto-allocation of the functions above

    xt::xtensor<double, N + 2> Strain() const;
    xt::xtensor<double, N + 2> Stress() const;
    xt::xtensor<double, N + 4> Tangent() const;

    // Get copy or reference to the underlying model at on point

    auto getElastic(const std::array<size_t, N>& index) const;
    auto* refElastic(const std::array<size_t, N>& index);

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
};

} // namespace Cartesian3d
} // namespace GMatElastic

#include "Cartesian3d.hpp"
#include "Cartesian3d_Array.hpp"
#include "Cartesian3d_Elastic.hpp"

#endif
