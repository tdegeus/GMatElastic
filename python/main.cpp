/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

// Enable basic assertions on matrix shape
// (doesn't cost a lot of time, but avoids segmentation faults)
#define GMATELASTIC_ENABLE_ASSERT

#include <GMatElastic/Cartesian3d.h>

namespace py = pybind11;

template <class S, class T>
auto construct_Array(T& self)
{
    namespace SM = GMatElastic::Cartesian3d;

    self.def(py::init<std::array<size_t, S::rank>>(), "Array of material points.", py::arg("shape"))
        .def(py::init<std::array<size_t, S::rank>, double, double>(),
             "Array of material points.",
             py::arg("shape"),
             py::arg("K"),
             py::arg("G"))

        .def("shape", &S::shape, "Shape of array.")
        .def("I2", &S::I2, "Array with 2nd-order unit tensors.")
        .def("II", &S::II, "Array with 4th-order tensors = dyadic(I2, I2).")
        .def("I4", &S::I4, "Array with 4th-order unit tensors.")
        .def("I4rt", &S::I4rt, "Array with 4th-order right-transposed unit tensors.")
        .def("I4s", &S::I4s, "Array with 4th-order symmetric projection tensors.")
        .def("I4d", &S::I4d, "Array with 4th-order deviatoric projection tensors.")
        .def("K", &S::K, "Array with bulk moduli.")
        .def("G", &S::G, "Array with shear moduli.")
        .def("type", &S::type, "Array with material types.")
        .def("isElastic", &S::isElastic, "Boolean-matrix: true for Elastic.")

        .def(
            "setElastic",
            &S::setElastic,
            "Set specific entries 'Elastic'.",
            py::arg("I"),
            py::arg("K"),
            py::arg("G"))

        .def("setStrain", &S::setStrain, "Set strain tensors.", py::arg("Eps"))
        .def("Strain", &S::Strain, "Get strain tensors.")
        .def("Stress", &S::Stress, "Get stress tensors.")
        .def("Tangent", &S::Tangent, "Get stiffness tensors.")
        .def("getElastic", &S::getElastic, "Returns underlying Elastic model.")

        .def("__repr__", [](const S&) { return "<GMatElastic.Cartesian3d.Array>"; });
}

template <class S, class T>
void add_deviatoric_overloads(T& module)
{
    module.def(
        "Deviatoric",
        static_cast<S (*)(const S&)>(&GMatElastic::Cartesian3d::Deviatoric<S>),
        "Deviatoric part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_hydrostatic_overloads(T& module)
{
    module.def(
        "Hydrostatic",
        static_cast<R (*)(const S&)>(&GMatElastic::Cartesian3d::Hydrostatic<S>),
        "Hydrostatic part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_epseq_overloads(T& module)
{
    module.def(
        "Epseq",
        static_cast<R (*)(const S&)>(
            &GMatElastic::Cartesian3d::Epseq<S>),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_sigeq_overloads(T& module)
{
    module.def(
        "Sigeq",
        static_cast<R (*)(const S&)>(
            &GMatElastic::Cartesian3d::Sigeq<S>),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"));
}

PYBIND11_MODULE(GMatElastic, m)
{

    m.doc() = "Linear elastic material model";

    // -----------------------
    // GMatElastic.Cartesian3d
    // -----------------------

    py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

    namespace SM = GMatElastic::Cartesian3d;

    // Unit tensors

    sm.def("I2", &SM::I2, "Second order unit tensor.");
    sm.def("II", &SM::II, "Fourth order tensor with the result of the dyadic product II.");
    sm.def("I4", &SM::I4, "Fourth order unit tensor.");
    sm.def("I4rt", &SM::I4rt, "Fourth right-transposed order unit tensor.");
    sm.def("I4s", &SM::I4s, "Fourth order symmetric projection tensor.");
    sm.def("I4d", &SM::I4d, "Fourth order deviatoric projection tensor.");

    // Tensor algebra

    add_deviatoric_overloads<xt::xtensor<double, 4>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 3>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 2>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);
    add_epseq_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_epseq_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_epseq_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);

    // Material point: Elastic

    py::class_<SM::Elastic>(sm, "Elastic")

        .def(py::init<double, double>(), "Linear elastic material point.", py::arg("K"), py::arg("G"))

        .def("K", &SM::Elastic::K, "Returns the bulk modulus.")
        .def("G", &SM::Elastic::G, "Returns the shear modulus.")
        .def("setStrain", &SM::Elastic::setStrain<xt::xtensor<double, 2>>, "Set strain tensor.")
        .def("Strain", &SM::Elastic::Strain, "Returns strain tensor.")
        .def("Stress", &SM::Elastic::Stress, "Returns stress tensor.")
        .def("Tangent", &SM::Elastic::Tangent, "Returns tangent stiffness.")

        .def("__repr__", [](const SM::Elastic&) { return "<GMatElastic.Cartesian3d.Elastic>"; });

    // Material identifier

    py::module smm = sm.def_submodule("Type", "Type enumerator");

    py::enum_<SM::Type::Value>(smm, "Type")
        .value("Unset", SM::Type::Unset)
        .value("Elastic", SM::Type::Elastic)
        .export_values();

    // Array

    py::class_<SM::Array<1>> array1d(sm, "Array1d");
    py::class_<SM::Array<2>> array2d(sm, "Array2d");
    py::class_<SM::Array<3>> array3d(sm, "Array3d");

    construct_Array<SM::Array<1>>(array1d);
    construct_Array<SM::Array<2>>(array2d);
    construct_Array<SM::Array<3>>(array3d);
}
