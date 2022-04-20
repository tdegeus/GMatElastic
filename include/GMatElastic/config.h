/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTIC_CONFIG_H
#define GMATELASTIC_CONFIG_H

#include <GMatTensor/config.h>

/**
All assertions are implementation as:

    GMATELASTIC_ASSERT(...)

They can be enabled by:

    #define GMATELASTIC_ENABLE_ASSERT

(before including GMatElastic).
The advantage is that:

-   File and line-number are displayed if the assertion fails.
-   Assertions can be enabled/disabled independently from those of other libraries.

\throw std::runtime_error
*/
#ifdef GMATELASTIC_ENABLE_ASSERT
#define GMATELASTIC_ASSERT(expr) GMATTENSOR_ASSERT_IMPL(expr, __FILE__, __LINE__)
#else
#define GMATELASTIC_ASSERT(expr)
#endif

/**
Linear elastic material model.
*/
namespace GMatElastic {

/**
Container type.
*/
namespace array_type {

#ifdef GMATELASTIC_USE_XTENSOR_PYTHON

/**
Fixed (static) rank array.
*/
template <typename T, size_t N>
using tensor = xt::pytensor<T, N>;

#else

/**
Fixed (static) rank array.
*/
template <typename T, size_t N>
using tensor = xt::xtensor<T, N>;

#endif

} // namespace array_type

} // namespace GMatElastic

#endif
