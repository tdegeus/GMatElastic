/**
\file
\copyright Copyright 2018. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pytensor.hpp>

#include <GMatElastic/Cartesian3d.h>
#include <GMatElastic/version.h>

namespace py = pybind11;

/**
Overrides the `__name__` of a module.
Classes defined by pybind11 use the `__name__` of the module as of the time they are defined,
which affects the `__repr__` of the class type objects.
*/
class ScopedModuleNameOverride {
public:
    explicit ScopedModuleNameOverride(py::module m, std::string name) : module_(std::move(m))
    {
        original_name_ = module_.attr("__name__");
        module_.attr("__name__") = name;
    }
    ~ScopedModuleNameOverride()
    {
        module_.attr("__name__") = original_name_;
    }

private:
    py::module module_;
    py::object original_name_;
};

template <class S, class T>
auto construct_Array(T& cls)
{
    cls.def(py::init<std::array<size_t, S::rank>>(), "Array of material points.", py::arg("shape"));
    cls.def(
        py::init<std::array<size_t, S::rank>, double, double>(),
        "Array of material points.",
        py::arg("shape"),
        py::arg("K"),
        py::arg("G"));

    cls.def("shape", &S::shape, "Shape of array.");
    cls.def("I2", &S::I2, "Array with 2nd-order unit tensors.");
    cls.def("II", &S::II, "Array with 4th-order tensors = dyadic(I2, I2).");
    cls.def("I4", &S::I4, "Array with 4th-order unit tensors.");
    cls.def("I4rt", &S::I4rt, "Array with 4th-order right-transposed unit tensors.");
    cls.def("I4s", &S::I4s, "Array with 4th-order symmetric projection tensors.");
    cls.def("I4d", &S::I4d, "Array with 4th-order deviatoric projection tensors.");
    cls.def("K", &S::K, "Array with bulk moduli.");
    cls.def("G", &S::G, "Array with shear moduli.");
    cls.def("type", &S::type, "Array with material types.");
    cls.def("isElastic", &S::isElastic, "Boolean-matrix: true for Elastic.");

    cls.def(
        "setElastic",
        static_cast<void (S::*)(const xt::pytensor<bool, S::rank>&, double, double)>(
            &S::template setElastic),
        "Set specific entries 'Elastic'.",
        py::arg("I"),
        py::arg("K"),
        py::arg("G"));

    cls.def(
        "setStrain",
        &S::template setStrain<xt::pytensor<double, S::rank + 2>>,
        "Set strain tensors.",
        py::arg("Eps"));

    cls.def("Strain", &S::Strain, "Get strain tensors.");

    cls.def(
        "strain",
        &S::template strain<xt::pytensor<double, S::rank + 2>>,
        "Get strain tensors for all points.");

    cls.def("Stress", &S::Stress, "Get stress tensors.");

    cls.def(
        "stress",
        &S::template stress<xt::pytensor<double, S::rank + 2>>,
        "Get stress tensors for all points.");

    cls.def("Tangent", &S::Tangent, "Get stiffness tensors.");

    cls.def(
        "tangent",
        &S::template tangent<xt::pytensor<double, S::rank + 4>>,
        "Get stiffness tensors.");

    cls.def(
        "refElastic",
        &S::refElastic,
        "Returns a reference to the underlying Elastic model.",
        py::return_value_policy::reference_internal);

    cls.def("__repr__", [](const S&) { return "<GMatElastic.Cartesian3d.Array>"; });
}

template <class R, class T, class M>
void add_Deviatoric(M& mod)
{
    mod.def(
        "Deviatoric",
        static_cast<R (*)(const T&)>(&GMatElastic::Cartesian3d::Deviatoric),
        "Deviatoric part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class T, class M>
void add_deviatoric(M& mod)
{
    mod.def(
        "deviatoric",
        static_cast<void (*)(const T&, R&)>(&GMatElastic::Cartesian3d::deviatoric),
        "Deviatoric part of a(n) (array of) tensor(s).",
        py::arg("A"),
        py::arg("ret"));
}

template <class R, class T, class M>
void add_Hydrostatic(M& mod)
{
    mod.def(
        "Hydrostatic",
        static_cast<R (*)(const T&)>(&GMatElastic::Cartesian3d::Hydrostatic),
        "Hydrostatic part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class T, class M>
void add_hydrostatic(M& mod)
{
    mod.def(
        "hydrostatic",
        static_cast<void (*)(const T&, R&)>(&GMatElastic::Cartesian3d::hydrostatic),
        "Hydrostatic part of a(n) (array of) tensor(s).",
        py::arg("A"),
        py::arg("ret"));
}

template <class R, class T, class M>
void add_Epseq(M& mod)
{
    mod.def(
        "Epseq",
        static_cast<R (*)(const T&)>(&GMatElastic::Cartesian3d::Epseq),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class T, class M>
void add_epseq(M& mod)
{
    mod.def(
        "epseq",
        static_cast<void (*)(const T&, R&)>(&GMatElastic::Cartesian3d::epseq),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"),
        py::arg("ret"));
}

template <class R, class T, class M>
void add_Sigeq(M& mod)
{
    mod.def(
        "Sigeq",
        static_cast<R (*)(const T&)>(&GMatElastic::Cartesian3d::Sigeq),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class T, class M>
void add_sigeq(M& mod)
{
    mod.def(
        "sigeq",
        static_cast<void (*)(const T&, R&)>(&GMatElastic::Cartesian3d::sigeq),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"),
        py::arg("ret"));
}

PYBIND11_MODULE(_GMatElastic, m)
{
    // Ensure members to display as `GMatElastic.X`
    // (not `GMatElastic._GMatElastic.X`)
    ScopedModuleNameOverride name_override(m, "GMatElastic");

    xt::import_numpy();

    m.doc() = "Linear elastic material model";

    m.def("version", &GMatElastic::version, "Return version string.");

    m.def("version_dependencies", &GMatElastic::version_dependencies, "Return list of strings.");

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

    add_Deviatoric<xt::pytensor<double, 4>, xt::pytensor<double, 4>>(sm);
    add_Deviatoric<xt::pytensor<double, 3>, xt::pytensor<double, 3>>(sm);
    add_Deviatoric<xt::pytensor<double, 2>, xt::pytensor<double, 2>>(sm);

    add_deviatoric<xt::pytensor<double, 4>, xt::pytensor<double, 4>>(sm);
    add_deviatoric<xt::pytensor<double, 3>, xt::pytensor<double, 3>>(sm);
    add_deviatoric<xt::pytensor<double, 2>, xt::pytensor<double, 2>>(sm);

    add_Hydrostatic<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    add_Hydrostatic<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    add_Hydrostatic<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    add_hydrostatic<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    add_hydrostatic<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    add_hydrostatic<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    add_Epseq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    add_Epseq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    add_Epseq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    add_epseq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    add_epseq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    add_epseq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    add_Sigeq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    add_Sigeq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    add_Sigeq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    add_sigeq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    add_sigeq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    add_sigeq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    // Material point: Elastic

    py::class_<SM::Elastic> cls(sm, "Elastic");

    cls.def(
        py::init<double, double>(), "Linear elastic material point.", py::arg("K"), py::arg("G"));

    cls.def("K", &SM::Elastic::K, "Returns the bulk modulus.");

    cls.def("G", &SM::Elastic::G, "Returns the shear modulus.");

    cls.def(
        "setStrain",
        &SM::Elastic::setStrain<xt::pytensor<double, 2>>,
        "Set current strain tensor.");

    cls.def("Strain", &SM::Elastic::Strain, "Returns strain tensor.");

    cls.def(
        "strain",
        &SM::Elastic::strain<xt::pytensor<double, 2>>,
        "Returns strain tensor.",
        py::arg("ret"));

    cls.def("Stress", &SM::Elastic::Stress, "Returns stress tensor.");

    cls.def(
        "stress",
        &SM::Elastic::stress<xt::pytensor<double, 2>>,
        "Returns stress tensor.",
        py::arg("ret"));

    cls.def("Tangent", &SM::Elastic::Tangent, "Returns tangent stiffness.");

    cls.def(
        "tangent",
        &SM::Elastic::tangent<xt::pytensor<double, 4>>,
        "Returns tangent stiffness.",
        py::arg("ret"));

    cls.def("energy", &SM::Elastic::energy, "Returns the energy, for last known strain.");

    cls.def("__repr__", [](const SM::Elastic&) { return "<GMatElastic.Cartesian3d.Elastic>"; });

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
