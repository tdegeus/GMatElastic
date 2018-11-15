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

PYBIND11_MODULE(GMatLinearElastic, m) {

m.doc() = "Linear elastic material model";

// ================================ GMatLinearElastic.Cartesian3d ==================================

{

// create sub-module
py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

// abbreviate name-space
namespace SM = GMatLinearElastic::Cartesian3d;

// -------------------------------------------------------------------------------------------------

sm.def("I"   , &SM::I   );
sm.def("II"  , &SM::II  );
sm.def("I4"  , &SM::I4  );
sm.def("I4rt", &SM::I4rt);
sm.def("I4s" , &SM::I4s );
sm.def("I4d" , &SM::I4d );

// -------------------------------------------------------------------------------------------------

py::class_<SM::Elastic>(sm, "Elastic")
  // constructor
  .def(
    py::init<double, double>(),
    "Linear elastic material point",
    py::arg("kappa"),
    py::arg("mu")
  )
  // methods
  .def("kappa"  , &SM::Elastic::kappa)
  .def("mu"     , &SM::Elastic::mu)
  .def("Sig"    , &SM::Elastic::Sig)
  .def("Tangent", &SM::Elastic::Tangent)
  // print to screen
  .def("__repr__", [](const SM::Elastic &){
    return "<GMatLinearElastic.Cartesian3d.Elastic>"; });

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
  .def("I"      , &SM::Matrix::I)
  .def("II"     , &SM::Matrix::II)
  .def("I4"     , &SM::Matrix::I4)
  .def("I4rt"   , &SM::Matrix::I4rt)
  .def("I4s"    , &SM::Matrix::I4s)
  .def("I4d"    , &SM::Matrix::I4d)
  .def("Sig"    , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Sig    , py::const_), py::arg("Eps"))
  .def("Tangent", py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Tangent, py::const_), py::arg("Eps"))
  // print to screen
  .def("__repr__", [](const SM::Matrix &){
    return "<GMatLinearElastic.Cartesian3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}

