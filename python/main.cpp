/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatLinearElastic

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

// Enable basic assertions on matrix shape
// (doesn't cost a lot of time, but avoids segmentation faults)
#define GMATLINEARELASTIC_ENABLE_ASSERT

// include library
#include "../include/GMatLinearElastic/Cartesian3d.h"

// abbreviate name-space
namespace py = pybind11;

// -------------------------------------------------------------------------------------------------

PYBIND11_MODULE(GMatLinearElastic, m) {

m.doc() = "Linear elastic material model";

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
    py::arg("K"),
    py::arg("G")
  )
  // methods
  .def("K"      , &SM::Elastic::K)
  .def("G"      , &SM::Elastic::G)
  .def("Stress" , &SM::Elastic::Stress)
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
    py::arg("K"),
    py::arg("G")
  )
  // methods
  .def("nelem"  , &SM::Matrix::nelem)
  .def("nip"    , &SM::Matrix::nip)
  .def("K"      , &SM::Matrix::K)
  .def("G"      , &SM::Matrix::G)
  .def("check"  , &SM::Matrix::check)
  .def("set"    , &SM::Matrix::set)
  .def("I"      , &SM::Matrix::I)
  .def("II"     , &SM::Matrix::II)
  .def("I4"     , &SM::Matrix::I4)
  .def("I4rt"   , &SM::Matrix::I4rt)
  .def("I4s"    , &SM::Matrix::I4s)
  .def("I4d"    , &SM::Matrix::I4d)
  .def("Stress" , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Stress , py::const_), py::arg("Eps"))
  .def("Tangent", py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Tangent, py::const_), py::arg("Eps"))
  // print to screen
  .def("__repr__", [](const SM::Matrix &){
    return "<GMatLinearElastic.Cartesian3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}


