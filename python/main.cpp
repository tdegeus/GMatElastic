/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatElastic

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

// Enable basic assertions on matrix shape
// (doesn't cost a lot of time, but avoids segmentation faults)
#define GMATELASTIC_ENABLE_ASSERT

// include library
#include "../include/GMatElastic/Cartesian3d.h"

// abbreviate name-space
namespace py = pybind11;

// -------------------------------------------------------------------------------------------------

PYBIND11_MODULE(GMatElastic, m) {

m.doc() = "Linear elastic material model";

// create submodule
py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

// abbreviate name-space
namespace SM = GMatElastic::Cartesian3d;

// -------------------------------------------------------------------------------------------------

sm.def("I", &SM::I);
sm.def("II", &SM::II);
sm.def("I4", &SM::I4);
sm.def("I4rt", &SM::I4rt);
sm.def("I4s", &SM::I4s);
sm.def("I4d", &SM::I4d);

// -------------------------------------------------------------------------------------------------

sm.def("Hydrostatic",
  py::overload_cast<const SM::Tensor2&>(&SM::Hydrostatic),
  "Hydrostatic part of a 2nd-order tensor",
  py::arg("A"));

sm.def("Deviatoric",
  py::overload_cast<const SM::Tensor2&>(&SM::Deviatoric),
  "Deviatoric",
  py::arg("A"));

sm.def("Epseq",
  py::overload_cast<const SM::Tensor2&>(&SM::Epseq),
  "Equivalent strain deviator",
  py::arg("Eps"));

sm.def("Sigeq",
  py::overload_cast<const SM::Tensor2&>(&SM::Sigeq),
  "Equivalent stress deviator",
  py::arg("Sig"));

// -------------------------------------------------------------------------------------------------

sm.def("Hydrostatic",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Hydrostatic),
  "Hydrostatic part of a 2nd-order tensor",
  py::arg("A"));

sm.def("Deviatoric",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Deviatoric),
  "Deviatoric",
  py::arg("A"));

sm.def("Epseq",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Epseq),
  "Equivalent strain deviator",
  py::arg("Eps"));

sm.def("Sigeq",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Sigeq),
  "Equivalent stress deviator",
  py::arg("Sig"));

// -------------------------------------------------------------------------------------------------

py::class_<SM::Elastic>(sm, "Elastic")

  .def(
    py::init<double, double>(),
    "Linear elastic material point",
    py::arg("K"),
    py::arg("G")
  )

  .def("K", &SM::Elastic::K)
  .def("G", &SM::Elastic::G)
  .def("Stress", &SM::Elastic::Stress, py::arg("Eps"))
  .def("Tangent", &SM::Elastic::Tangent, py::arg("Eps"))

  .def("__repr__", [](const SM::Elastic &){
    return "<GMatElastic.Cartesian3d.Elastic>"; });

// -------------------------------------------------------------------------------------------------

py::module smm = sm.def_submodule("Type", "Type enumerator");

py::enum_<SM::Type::Value>(smm, "Type")
    .value("Unset", SM::Type::Unset)
    .value("Elastic", SM::Type::Elastic)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::class_<SM::Matrix>(sm, "Matrix")

  .def(
    py::init<size_t, size_t>(),
    "Matrix of linear elastic material points",
    py::arg("nelem"),
    py::arg("nip")
  )

  .def(
    py::init<size_t, size_t, double, double>(),
    "Matrix of linear elastic material points",
    py::arg("nelem"),
    py::arg("nip"),
    py::arg("K"),
    py::arg("G"))

  .def("ndim", &SM::Matrix::ndim)
  .def("nelem", &SM::Matrix::nelem)
  .def("nip", &SM::Matrix::nip)

  .def("K", &SM::Matrix::K)
  .def("G", &SM::Matrix::G)

  .def("I", &SM::Matrix::I)
  .def("II", &SM::Matrix::II)
  .def("I4", &SM::Matrix::I4)
  .def("I4rt", &SM::Matrix::I4rt)
  .def("I4s", &SM::Matrix::I4s)
  .def("I4d", &SM::Matrix::I4d)

  .def("check", &SM::Matrix::check)

  .def("setElastic",
    &SM::Matrix::setElastic,
    py::arg("I"),
    py::arg("K"),
    py::arg("G"))

  .def("Stress", &SM::Matrix::Stress, py::arg("Eps"))
  .def("Tangent", &SM::Matrix::Tangent, py::arg("Eps"))

  .def("__repr__", [](const SM::Matrix &){
    return "<GMatElastic.Cartesian3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

