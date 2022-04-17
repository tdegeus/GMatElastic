/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTIC_CARTESIAN3D_H
#define GMATELASTIC_CARTESIAN3D_H

#include <GMatTensor/Cartesian3d.h>

#include "config.h"

namespace GMatElastic {

/**
Implementation in a 3-d Cartesian coordinate frame.

Note that for convenience this namespace include aliases to:
-   GMatTensor::Cartesian3d::Deviatoric()
-   GMatTensor::Cartesian3d::deviatoric()
-   GMatTensor::Cartesian3d::Hydrostatic()
-   GMatTensor::Cartesian3d::hydrostatic()
-   GMatTensor::Cartesian3d::I2()
-   GMatTensor::Cartesian3d::I4()
-   GMatTensor::Cartesian3d::I4d()
-   GMatTensor::Cartesian3d::I4rt()
-   GMatTensor::Cartesian3d::I4s()
-   GMatTensor::Cartesian3d::II()
-   GMatTensor::Cartesian3d::O2()
-   GMatTensor::Cartesian3d::O4()
*/
namespace Cartesian3d {

using GMatTensor::Cartesian3d::deviatoric;
using GMatTensor::Cartesian3d::Deviatoric;
using GMatTensor::Cartesian3d::hydrostatic;
using GMatTensor::Cartesian3d::Hydrostatic;
using GMatTensor::Cartesian3d::I2;
using GMatTensor::Cartesian3d::I4;
using GMatTensor::Cartesian3d::I4d;
using GMatTensor::Cartesian3d::I4rt;
using GMatTensor::Cartesian3d::I4s;
using GMatTensor::Cartesian3d::II;
using GMatTensor::Cartesian3d::O2;
using GMatTensor::Cartesian3d::O4;

/**
Von Mises equivalent strain: norm of strain deviator

\f$ \sqrt{\frac{2}{3} (dev(A))_{ij} (dev(A))_{ji}} \f$

To write to allocated data use epseq().

\param A [..., 3, 3] array.
\return [...] array.
*/
template <class T>
inline auto Epseq(const T& A) ->
    typename GMatTensor::detail::allocate<xt::get_rank<T>::value - 2, T>::type;

/**
Same as epseq(), but writes to externally allocated output.

\param A [..., 3, 3] array.
\param ret output [...] array
*/
template <class T, class U>
inline void epseq(const T& A, U& ret);

/**
Von Mises equivalent stress: norm of strain deviator

\f$ \sqrt{\frac{3}{2} (dev(A))_{ij} (dev(A))_{ji}} \f$

To write to allocated data use sigeq().

\param A [..., 3, 3] array.
\return [...] array.
*/
template <class T>
inline auto Sigeq(const T& A) ->
    typename GMatTensor::detail::allocate<xt::get_rank<T>::value - 2, T>::type;

/**
Same as Sigeq(), but writes to externally allocated output.

\param A [..., 3, 3] array.
\param ret output [...] array
*/
template <class T, class U>
inline void sigeq(const T& A, U& ret);

/**
Elastic material point.
*/
class Elastic {
public:
    Elastic() = default;

    /**
    Constructor.

    \param K Bulk modulus.
    \param G Shear modulus.
    */
    Elastic(double K, double G);

    /**
    \return Bulk modulus.
    */
    double K() const;

    /**
    \return Shear modulus.
    */
    double G() const;

    /**
    \return Current potential energy.
    */
    double energy() const;

    /**
    Set the current strain tensor.

    \param arg xtensor array [3, 3].
    */
    template <class T>
    void setStrain(const T& arg);

    /**
    Same as setStrain(), but reads from a pointer assuming row-major storage (no bound check).

    \param arg Pointer to array (xx, xy, xz, yx, yy, yz, zx, zy, zz).
    */
    template <class T>
    void setStrainPtr(const T* arg);

    /**
    Get the current strain tensor.

    \return [3, 3] array.
    */
    xt::xtensor<double, 2> Strain() const;

    /**
    Same as Strain(), but write to allocated data.

    \param ret xtensor array [3, 3], overwritten.
    */
    template <class T>
    void strain(T& ret) const;

    /**
    Same as Strain(), but write to a pointer assuming row-major storage (no bound check).

    \param ret Pointer to array (xx, xy, xz, yx, yy, yz, zx, zy, zz), overwritten.
    */
    template <class T>
    void strainPtr(T* ret) const;

    /**
    Get the current stress tensor.

    \return [3, 3] array.
    */
    xt::xtensor<double, 2> Stress() const;

    /**
    Same as Stress(), but write to allocated data.

    \param ret xtensor array [3, 3], overwritten.
    */
    template <class T>
    void stress(T& ret) const;

    /**
    Same as Stress(), but write to a pointer assuming row-major storage (no bound check).

    \param ret Pointer to array (xx, xy, xz, yx, yy, yz, zx, zy, zz), overwritten.
    */
    template <class T>
    void stressPtr(T* ret) const;

    /**
    Get the tangent tensor (strain independent).

    \return [3, 3, 3, 3] array.
    */
    xt::xtensor<double, 4> Tangent() const;

    /**
    Same as Tangent(), but write to allocated data.

    \param ret xtensor array [3, 3, 3, 3], overwritten.
    */
    template <class T>
    void tangent(T& ret) const;

    /**
    Same as Tangent(), but write to a pointer assuming row-major storage (no bound check).

    \param ret Pointer to array of size 3 * 3 * 3 * 3, overwritten.
    */
    template <class T>
    void tangentPtr(T* ret) const;

private:
    double m_K; ///< bulk modulus
    double m_G; ///< shear modulus
    std::array<double, 9> m_Eps; ///< strain tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    std::array<double, 9> m_Sig; ///< stress tensor ,,
};

/**
Material identifier.
*/
struct Type {
    /**
    Type value.
    */
    enum Value {
        Unset, ///< Unset
        Elastic, ///< See Elastic
    };
};

/**
Array of material points.

\tparam N Rank of the array.
*/
template <size_t N>
class Array : public GMatTensor::Cartesian3d::Array<N> {
public:
    using GMatTensor::Cartesian3d::Array<N>::rank;

    Array() = default;

    /**
    Basic constructor.
    Note that before usage material properties still have to be assigned to all items.
    This can be done per item or by groups of items, using:
    -   setElastic()

    \param shape The shape of the array.
    */
    Array(const std::array<size_t, N>& shape);

    /**
    Construct homogeneous system.

    \param shape The shape of the array.
    \param K Bulk modulus.
    \param G Shear modulus.
    */
    Array(const std::array<size_t, N>& shape, double K, double G);

    /**
    \return Type-id per item. Follows order set in Type.
    */
    xt::xtensor<size_t, N> type() const;

    /**
    \return Per item, 1 if Elastic, otherwise 0.
    */
    xt::xtensor<bool, N> isElastic() const;

    /**
    \return Bulk modulus per item.
    */
    xt::xtensor<double, N> K() const;

    /**
    \return Shear modulus per item.
    */
    xt::xtensor<double, N> G() const;

    /**
    Set all items Elastic, specifying material parameters per item.

    \tparam L e.g. `xt::xtensor<bool, N>`
    \param I Per item, ``true`` to set Elastic, ``false`` to skip.
    \param K Bulk modulus.
    \param G Shear modulus.
    */
    template <class L>
    void setElastic(const L& I, double K, double G);

    /**
    Set strain tensors.

    \tparam T e.g. `xt::xtensor<double, N + 2>`
    \param arg Strain tensor per item [shape(), 3, 3].
    */
    template <class T>
    void setStrain(const T& arg);

    /**
    \return Strain tensor per item [shape(), 3, 3].
    */
    xt::xtensor<double, N + 2> Strain() const;

    /**
    Same as Strain(), but write to allocated data.

    \tparam R e.g. `xt::xtensor<double, N + 2>`
    \param ret [shape(), 3, 3], overwritten.
    */
    template <class R>
    void strain(R& ret) const;

    /**
    \return Stress tensor per item [shape(), 3, 3].
    */
    xt::xtensor<double, N + 2> Stress() const;

    /**
    Same as Stress(), but write to allocated data.

    \tparam R e.g. `xt::xtensor<double, N + 2>`
    \param ret [shape(), 3, 3], overwritten.
    */
    template <class R>
    void stress(R& ret) const;

    /**
    \return Tangent tensor per item [shape(), 3, 3, 3, 3].
    */
    xt::xtensor<double, N + 4> Tangent() const;

    /**
    Same as Tangent(), but write to allocated data.

    \tparam R e.g. `xt::xtensor<double, N + 4>`
    \param ret [shape(), 3, 3, 3, 3], overwritten.
    */
    template <class R>
    void tangent(R& ret) const;

    /**
    Reference to the underlying Elastic model of an item.

    \param index The index of the item.
    \return Reference to the model.
    */
    Elastic& refElastic(const std::array<size_t, N>& index);

    /**
    Constant reference to the underlying Elastic model of an item.

    \param index The index of the item.
    \return Reference to the model.
    */
    const Elastic& crefElastic(const std::array<size_t, N>& index) const;

private:
    /**
    Elastic material vectors: each item has one entry in one of the material vectors.
    */
    std::vector<Elastic> m_Elastic;

    /**
    Type of each entry, see Type.
    */
    xt::xtensor<size_t, N> m_type;

    /**
    Index in the relevant material vector (#m_Elastic)
    */
    xt::xtensor<size_t, N> m_index;

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
