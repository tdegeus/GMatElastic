/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_CARTESIAN3D_H
#define GMATELASTIC_CARTESIAN3D_H

#include "config.h"

namespace GMatElastic {
namespace Cartesian3d {

// Alias

#if defined(_WIN32) || defined(_WIN64)
    using Tensor2 = xt::xtensor<double, 2>;
    using Tensor4 = xt::xtensor<double, 4>;
#else
    using Tensor2 = xt::xtensor_fixed<double, xt::xshape<3, 3>>;
    using Tensor4 = xt::xtensor_fixed<double, xt::xshape<3, 3, 3, 3>>;
#endif

// Unit tensors

inline Tensor2 I2();
inline Tensor4 II();
inline Tensor4 I4();
inline Tensor4 I4rt();
inline Tensor4 I4s();
inline Tensor4 I4d();

// Tensor decomposition

template <class T, class U>
inline void hydrostatic(const T& A, U& B);

template <class T>
inline auto Hydrostatic(const T& A);

template <class T, class U>
inline void deviatoric(const T& A, U& B);

template <class T>
inline auto Deviatoric(const T& A);

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
    Tensor2 Stress() const;
    Tensor4 Tangent() const;

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

template <size_t rank>
class Array
{
public:
    // Constructors

    Array() = default;
    Array(const std::array<size_t, rank>& shape);
    Array(const std::array<size_t, rank>& shape, double K, double G);

    // Shape

    std::array<size_t, rank> shape() const;

    // Type

    xt::xtensor<size_t, rank> type() const;
    xt::xtensor<size_t, rank> isElastic() const;

    // Parameters

    xt::xtensor<double, rank> K() const;
    xt::xtensor<double, rank> G() const;

    // Array of unit tensors

    xt::xtensor<double, rank + 2> I2() const;
    xt::xtensor<double, rank + 4> II() const;
    xt::xtensor<double, rank + 4> I4() const;
    xt::xtensor<double, rank + 4> I4rt() const;
    xt::xtensor<double, rank + 4> I4s() const;
    xt::xtensor<double, rank + 4> I4d() const;

    // Check that a type has been set everywhere (throws if unset points are found)

    void check() const;

    // Set parameters for a batch of points

    void setElastic(const xt::xtensor<size_t, rank>& I, double K, double G);

    // Set strain tensor, get the response

    void setStrain(const xt::xtensor<double, rank + 2>& Eps);
    void stress(xt::xtensor<double, rank + 2>& Sig) const;
    void tangent(xt::xtensor<double, rank + 4>& C) const;

    // Auto-allocation of the functions above

    xt::xtensor<double, rank + 2> Stress() const;
    xt::xtensor<double, rank + 4> Tangent() const;

private:
    // Material vectors
    std::vector<Elastic> m_Elastic;

    // Identifiers for each matrix entry
    xt::xtensor<size_t, rank> m_type;  // type (e.g. "Type::Elastic")
    xt::xtensor<size_t, rank> m_index; // index from the relevant material vector (e.g. "m_Elastic")

    // Shape
    static const size_t m_ndim = 3;
    size_t m_size;
    std::array<size_t, rank> m_shape;
    std::array<size_t, rank + 2> m_shape_tensor2;
    std::array<size_t, rank + 4> m_shape_tensor4;

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
