/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatElastic

================================================================================================= */

#ifndef GMATELASTIC_H
#define GMATELASTIC_H

// -------------------------------------------------------------------------------------------------

// use "M_PI" from "math.h"
#define _USE_MATH_DEFINES

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>

// -------------------------------------------------------------------------------------------------

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// -------------------------------------------------------------------------------------------------

#ifndef NDEBUG
#define GMATELASTIC_ENABLE_ASSERT
#endif

#ifdef GMATELASTIC_ENABLE_ASSERT
#define GMATELASTIC_ASSERT(expr) GMATELASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)
#define GMATELASTIC_ASSERT_IMPL(expr, file, line)                                                                         \
    if (!(expr))                                                                                                          \
    {                                                                                                                     \
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": assertion failed (" #expr ") \n\t"); \
    }
#else
#define GMATELASTIC_ASSERT(expr)
#endif

// -------------------------------------------------------------------------------------------------

#define GMATELASTIC_WORLD_VERSION 0
#define GMATELASTIC_MAJOR_VERSION 2
#define GMATELASTIC_MINOR_VERSION 1

#define GMATELASTIC_VERSION_AT_LEAST(x,y,z) \
  (GMATELASTIC_WORLD_VERSION>x || (GMATELASTIC_WORLD_VERSION>=x && \
  (GMATELASTIC_MAJOR_VERSION>y || (GMATELASTIC_MAJOR_VERSION>=y && \
                                         GMATELASTIC_MINOR_VERSION>=z))))

#define GMATELASTIC_VERSION(x,y,z) \
  (GMATELASTIC_WORLD_VERSION==x && \
   GMATELASTIC_MAJOR_VERSION==y && \
   GMATELASTIC_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#endif
