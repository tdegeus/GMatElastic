/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef GMATLINEARELASTIC_H
#define GMATLINEARELASTIC_H

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

#define GMATLINEARELASTIC_WORLD_VERSION 0
#define GMATLINEARELASTIC_MAJOR_VERSION 1
#define GMATLINEARELASTIC_MINOR_VERSION 1

#define GMATLINEARELASTIC_VERSION_AT_LEAST(x,y,z) \
  (GMATLINEARELASTIC_WORLD_VERSION>x || (GMATLINEARELASTIC_WORLD_VERSION>=x && \
  (GMATLINEARELASTIC_MAJOR_VERSION>y || (GMATLINEARELASTIC_MAJOR_VERSION>=y && \
                                         GMATLINEARELASTIC_MINOR_VERSION>=z))))

#define GMATLINEARELASTIC_VERSION(x,y,z) \
  (GMATLINEARELASTIC_WORLD_VERSION==x && \
   GMATLINEARELASTIC_MAJOR_VERSION==y && \
   GMATLINEARELASTIC_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#endif
