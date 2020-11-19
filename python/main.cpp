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

template <size_t rank, class T>
auto add_common_members_array(T& self)
{
    namespace SM = GMatElastic::Cartesian3d;

    self.def(py::init<std::array<size_t, rank>>(), "Array of material points.", py::arg("shape"))

        .def("shape", &SM::Array<rank>::shape, "Shape of array.")

        .def("K", &SM::Array<rank>::K, "Array with bulk moduli.")

        .def("G", &SM::Array<rank>::G, "Array with shear moduli.")

        .def("I2", &SM::Array<rank>::I2, "Array with 2nd-order unit tensors.")

        .def("II", &SM::Array<rank>::II, "Array with 4th-order tensors = dyadic(I2, I2).")

        .def("I4", &SM::Array<rank>::I4, "Array with 4th-order unit tensors.")

        .def("I4rt", &SM::Array<rank>::I4rt, "Array with 4th-order right-transposed unit tensors.")

        .def("I4s", &SM::Array<rank>::I4s, "Array with 4th-order symmetric projection tensors.")

        .def("I4d", &SM::Array<rank>::I4d, "Array with 4th-order deviatoric projection tensors.")

        .def("type", &SM::Array<rank>::type, "Array with material types.")

        .def("isElastic", &SM::Array<rank>::isElastic, "Boolean-matrix: true for Elastic.")

        .def("check", &SM::Array<rank>::check, "Throws if any unset point is found.")

        .def(
            "setElastic",
            py::overload_cast<const xt::xtensor<size_t, rank>&, double, double>(
                &SM::Array<rank>::setElastic),
            "Set specific entries 'Elastic'.",
            py::arg("I"),
            py::arg("K"),
            py::arg("G"))

        .def(
            "setStrain",
            &SM::Array<rank>::setStrain,
            "Set matrix of strain tensors.",
            py::arg("Eps"))

        .def(
            "Stress",
            &SM::Array<rank>::Stress,
            "Returns matrix of stress tensors, given the current strain.")

        .def(
            "Tangent",
            &SM::Array<rank>::Tangent,
            "Returns matrices of stress tangent stiffness tensors, given the current strain.")

        .def("__repr__", [](const SM::Array<rank>&) {
            return "<GMatElastic.Cartesian3d.Array>";
        });
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

    sm.def(
        "Deviatoric",
        static_cast<xt::xtensor<double, 4> (*)(const xt::xtensor<double, 4>&)>(
            &SM::Deviatoric<xt::xtensor<double, 4>>),
        "Deviatoric part of a 2nd-order tensor. Returns matrix of 2nd-order tensors.",
        py::arg("A"));

    sm.def(
        "Deviatoric",
        static_cast<xt::xtensor<double, 3> (*)(const xt::xtensor<double, 3>&)>(
            &SM::Deviatoric<xt::xtensor<double, 3>>),
        "Deviatoric part of a 2nd-order tensor. Returns list of 2nd-order tensors.",
        py::arg("A"));

    sm.def(
        "Deviatoric",
        static_cast<xt::xtensor<double, 2> (*)(const xt::xtensor<double, 2>&)>(
            &SM::Deviatoric<xt::xtensor<double, 2>>),
        "Deviatoric part of a 2nd-order tensor. Returns 2nd-order tensor.",
        py::arg("A"));

    sm.def(
        "Hydrostatic",
        static_cast<xt::xtensor<double, 2> (*)(const xt::xtensor<double, 4>&)>(
            &SM::Hydrostatic<xt::xtensor<double, 4>>),
        "Hydrostatic part of a 2nd-order tensor. Returns matrix (of scalars).",
        py::arg("A"));

    sm.def(
        "Hydrostatic",
        static_cast<xt::xtensor<double, 1> (*)(const xt::xtensor<double, 3>&)>(
            &SM::Hydrostatic<xt::xtensor<double, 3>>),
        "Hydrostatic part of a 2nd-order tensor. Returns list (of scalars).",
        py::arg("A"));

    sm.def(
        "Hydrostatic",
        static_cast<xt::xtensor<double, 0> (*)(const xt::xtensor<double, 2>&)>(
            &SM::Hydrostatic<xt::xtensor<double, 2>>),
        "Hydrostatic part of a 2nd-order tensor. Returns scalar.",
        py::arg("A"));

    sm.def(
        "Epseq",
        static_cast<xt::xtensor<double, 2> (*)(const xt::xtensor<double, 4>&)>(
            &SM::Epseq<xt::xtensor<double, 4>>),
        "Equivalent strain. Returns matrix (of scalars).",
        py::arg("A"));

    sm.def(
        "Epseq",
        static_cast<xt::xtensor<double, 1> (*)(const xt::xtensor<double, 3>&)>(
            &SM::Epseq<xt::xtensor<double, 3>>),
        "Equivalent strain. Returns list (of scalars).",
        py::arg("A"));

    sm.def(
        "Epseq",
        static_cast<xt::xtensor<double, 0> (*)(const xt::xtensor<double, 2>&)>(
            &SM::Epseq<xt::xtensor<double, 2>>),
        "Equivalent strain. Returns scalar.",
        py::arg("A"));

    sm.def(
        "Sigeq",
        static_cast<xt::xtensor<double, 2> (*)(const xt::xtensor<double, 4>&)>(
            &SM::Sigeq<xt::xtensor<double, 4>>),
        "Equivalent stress. Returns matrix (of scalars).",
        py::arg("A"));

    sm.def(
        "Sigeq",
        static_cast<xt::xtensor<double, 1> (*)(const xt::xtensor<double, 3>&)>(
            &SM::Sigeq<xt::xtensor<double, 3>>),
        "Equivalent stress. Returns list (of scalars).",
        py::arg("A"));

    sm.def(
        "Sigeq",
        static_cast<xt::xtensor<double, 0> (*)(const xt::xtensor<double, 2>&)>(
            &SM::Sigeq<xt::xtensor<double, 2>>),
        "Equivalent stress. Returns scalar.",
        py::arg("A"));

    // Material point: Elastic

    py::class_<SM::Elastic>(sm, "Elastic")

        .def(py::init<double, double>(), "Linear elastic material point.", py::arg("K"), py::arg("G"))

        .def("K", &SM::Elastic::K, "Returns the bulk modulus.")

        .def("G", &SM::Elastic::G, "Returns the shear modulus.")

        .def("setStrain", &SM::Elastic::setStrain<SM::Tensor2>, "Set current strain tensor.")

        .def("Stress", &SM::Elastic::Stress, "Returns stress tensor, for last known strain.")

        .def(
            "Tangent",
            &SM::Elastic::Tangent,
            "Returns stress and tangent stiffness tensors, for last known strain.")

        .def("__repr__", [](const SM::Elastic&) { return "<GMatElastic.Cartesian3d.Elastic>"; });

    // Material identifier

    py::module smm = sm.def_submodule("Type", "Type enumerator");

    py::enum_<SM::Type::Value>(smm, "Type")
        .value("Unset", SM::Type::Unset)
        .value("Elastic", SM::Type::Elastic)
        .export_values();

    // Array

    py::class_<SM::Array<1>> array1d(sm, "Array1d");
    add_common_members_array<1>(array1d);

    py::class_<SM::Array<2>> array2d(sm, "Array2d");
    add_common_members_array<2>(array2d);

    py::class_<SM::Array<3>> array3d(sm, "Array3d");
    add_common_members_array<3>(array3d);

}
