#define CATCH_CONFIG_MAIN
#include <GMatElastic/Cartesian3d.h>
#include <GMatTensor/Cartesian3d.h>
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1e-12));

namespace GM = GMatElastic::Cartesian3d;
namespace GT = GMatTensor::Cartesian3d;

TEST_CASE("GMatElastic::Cartesian3d", "Cartesian3d.h")
{
    SECTION("Epseq - Tensor2")
    {
        auto A = GT::O2();
        A(0, 1) = 1.0;
        A(1, 0) = 1.0;
        REQUIRE(GM::Epseq(A)() == Approx(2.0 / std::sqrt(3.0)));
    }

    SECTION("Epseq - List")
    {
        auto A = GT::O2();
        A(0, 1) = 1.0;
        A(1, 0) = 1.0;
        auto M = xt::xtensor<double, 3>::from_shape({3, 3, 3});
        auto R = xt::xtensor<double, 1>::from_shape({M.shape(0)});
        for (size_t i = 0; i < M.shape(0); ++i) {
            xt::view(M, i, xt::all(), xt::all()) = static_cast<double>(i) * A;
            R(i) = static_cast<double>(i) * 2.0 / std::sqrt(3.0);
        }
        REQUIRE(xt::allclose(GM::Epseq(M), R));
    }

    SECTION("Epseq - Matrix")
    {
        auto A = GT::O2();
        A(0, 1) = 1.0;
        A(1, 0) = 1.0;
        auto M = xt::xtensor<double, 4>::from_shape({3, 4, 3, 3});
        auto R = xt::xtensor<double, 2>::from_shape({M.shape(0), M.shape(1)});
        for (size_t i = 0; i < M.shape(0); ++i) {
            for (size_t j = 0; j < M.shape(1); ++j) {
                xt::view(M, i, j, xt::all(), xt::all()) =
                    static_cast<double>(i * M.shape(1) + j) * A;
                R(i, j) = static_cast<double>(i * M.shape(1) + j) * 2.0 / std::sqrt(3.0);
            }
        }
        REQUIRE(xt::allclose(GM::Epseq(M), R));
    }

    SECTION("Sigeq - Tensor2")
    {
        auto A = GT::O2();
        A(0, 1) = 1.0;
        A(1, 0) = 1.0;
        REQUIRE(GM::Sigeq(A)() == Approx(std::sqrt(3.0)));
    }

    SECTION("Sigeq - List")
    {
        auto A = GT::O2();
        A(0, 1) = 1.0;
        A(1, 0) = 1.0;
        auto M = xt::xtensor<double, 3>::from_shape({3, 3, 3});
        auto R = xt::xtensor<double, 1>::from_shape({M.shape(0)});
        for (size_t i = 0; i < M.shape(0); ++i) {
            xt::view(M, i, xt::all(), xt::all()) = static_cast<double>(i) * A;
            R(i) = static_cast<double>(i) * std::sqrt(3.0);
        }
        REQUIRE(xt::allclose(GM::Sigeq(M), R));
    }

    SECTION("Sigeq - Matrix")
    {
        auto A = GT::O2();
        A(0, 1) = 1.0;
        A(1, 0) = 1.0;
        auto M = xt::xtensor<double, 4>::from_shape({3, 4, 3, 3});
        auto R = xt::xtensor<double, 2>::from_shape({M.shape(0), M.shape(1)});
        for (size_t i = 0; i < M.shape(0); ++i) {
            for (size_t j = 0; j < M.shape(1); ++j) {
                xt::view(M, i, j, xt::all(), xt::all()) =
                    static_cast<double>(i * M.shape(1) + j) * A;
                R(i, j) = static_cast<double>(i * M.shape(1) + j) * std::sqrt(3.0);
            }
        }
        REQUIRE(xt::allclose(GM::Sigeq(M), R));
    }

    SECTION("Elastic - stress")
    {
        double K = 12.3;
        double G = 45.6;
        double gamma = 0.02;
        double epsm = 0.12;

        xt::xtensor<double, 2> Eps = {{epsm, gamma, 0.0}, {gamma, epsm, 0.0}, {0.0, 0.0, epsm}};

        xt::xtensor<double, 2> Sig = {
            {3.0 * K * epsm, 2.0 * G * gamma, 0.0},
            {2.0 * G * gamma, 3.0 * K * epsm, 0.0},
            {0.0, 0.0, 3.0 * K * epsm}};

        GM::Elastic mat(K, G);
        mat.setStrain(Eps);

        REQUIRE(xt::allclose(mat.Stress(), Sig));
    }

    SECTION("Tangent")
    {
        double K = 12.3;
        double G = 45.6;

        xt::xtensor<double, 2> Eps = xt::random::randn<double>({3, 3});
        xt::xtensor<double, 4> Is = GM::I4s();
        Eps = GT::A4_ddot_B2(Is, Eps);

        GM::Elastic mat(K, G);
        mat.setStrain(Eps);
        auto Sig = mat.Stress();
        auto C = mat.Tangent();
        REQUIRE(xt::allclose(GT::A4_ddot_B2(C, Eps), Sig));
    }

    SECTION("Array")
    {
        double K = 12.3;
        double G = 45.6;
        double gamma = 0.02;
        double epsm = 0.12;

        xt::xtensor<double, 2> Eps = {{epsm, gamma, 0.0}, {gamma, epsm, 0.0}, {0.0, 0.0, epsm}};

        xt::xtensor<double, 2> Sig = {
            {3.0 * K * epsm, 2.0 * G * gamma, 0.0},
            {2.0 * G * gamma, 3.0 * K * epsm, 0.0},
            {0.0, 0.0, 3.0 * K * epsm}};

        size_t nelem = 3;
        size_t nip = 2;
        size_t ndim = 3;

        GM::Array<2> mat({nelem, nip});

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem, nip});
            mat.setElastic(I, K, G);
        }

        xt::xtensor<double, 4> eps = xt::empty<double>({nelem, nip, ndim, ndim});
        xt::xtensor<double, 4> sig = xt::empty<double>({nelem, nip, ndim, ndim});

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                double fac = static_cast<double>((e + 1) * nip + (q + 1));
                xt::view(eps, e, q) = fac * Eps;
                xt::view(sig, e, q) = fac * Sig;
            }
        }

        mat.setStrain(eps);

        REQUIRE(xt::allclose(mat.Stress(), sig));
    }

    SECTION("Array - Model")
    {
        double K = 12.3;
        double G = 45.6;
        double gamma = 0.02;
        double epsm = 0.12;

        xt::xtensor<double, 2> Eps = {{epsm, gamma, 0.0}, {gamma, epsm, 0.0}, {0.0, 0.0, epsm}};

        xt::xtensor<double, 2> Sig = {
            {3.0 * K * epsm, 2.0 * G * gamma, 0.0},
            {2.0 * G * gamma, 3.0 * K * epsm, 0.0},
            {0.0, 0.0, 3.0 * K * epsm}};

        size_t nelem = 3;
        size_t nip = 2;

        GM::Array<2> mat({nelem, nip});

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem, nip});
            mat.setElastic(I, K, G);
        }

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                double fac = static_cast<double>((e + 1) * nip + (q + 1));
                auto model = mat.getElastic({e, q});
                model.setStrain(xt::eval(fac * Eps));
                REQUIRE(xt::allclose(model.Stress(), fac * Sig));
            }
        }
    }
}
