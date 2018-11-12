/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

#include "Cartesian3d.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;

// ======================================= GMatLinearElastic =======================================

PYBIND11_MODULE(GMatLinearELlastic, m) {

m.doc() = "Linear elastic material model";

// ================================ GMatLinearELlastic::Cartesian3d ================================

{

// create submodule
py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

// abbreviate name-space
namespace SM = GMatLinearElastic::Cartesian3d;

// -------------------------------------------------------------------------------------------------

py::class_<SM::Matrix>(sm, "Matrix")
  // constructor
  .def(
    py::init<size_t, size_t>(),
    "Matrix of linear elastic material points",
    py::arg("nelem"),
    py::arg("nip")
  )
  // constructor
  .def(
    py::init<size_t, size_t, double, double>(),
    "Matrix of linear elastic material points",
    py::arg("nelem"),
    py::arg("nip"),
    py::arg("kappa"),
    py::arg("mu")
  )
  // methods
  .def("nelem"  , &SM::Matrix::nelem)
  .def("nip"    , &SM::Matrix::nip)
  .def("kappa"  , &SM::Matrix::kappa)
  .def("mu"     , &SM::Matrix::mu)
  .def("check"  , &SM::Matrix::check)
  .def("set"    , &SM::Matrix::set)
  .def("Sig"    , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Sig    , py::const_), py::arg("Eps"))
  .def("Tangent", py::overload_cast<>(&SM::Matrix::Tangent, py::const_))
  // print to screen
  .def("__repr__", [](const SM::Matrix &){
    return "GMat<LinearElastic.Cartesian3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}

