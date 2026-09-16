/**
 * @file well_natural_test.cpp
 * @brief 单元测试：验证 `well_natural` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <natural/physics/well_source.hpp>
#include <well/specification.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using Config = MPMC::CompositionalModelConfig<6, true, true>;
using Indices = MPMC::ScalarIndices<Config>;

void near(double value, double expected, double tolerance, const char *message)
{
    if (std::abs(value - expected) > tolerance)
        throw std::runtime_error(message);
}
}

int main()
{
    MPMC::WellSpecification<Indices, long long> spec;
    spec.type = MPMC::WellType::Producer;
    spec.control = MPMC::WellControl::TotalRate;
    spec.target = 12.0;
    spec.bhpCellId = 2;
    spec.perforations = {{2, 1.0e-12}, {7, 2.0e-12}};
    spec.validate(10);
    near(spec.signedTarget(), -12.0, 0.0,
         "Natural bridge producer target");

    std::array<double, Indices::numPhases> pressure{100.0, 100.0, 100.0};
    std::array<double, Indices::numPhases> density{800.0, 100.0, 1000.0};
    std::array<double, Indices::numPhases> mobility{0.1, 0.2, 0.3};
    std::array<double, Indices::numPhases> surfaceDensity{800.0, 2.0, 1000.0};
    std::array<double, Indices::numPhases> phaseFraction{1.0, 0.0, 0.0};
    std::array<double, Indices::numComponents> injectionComposition{1.0, 0, 0, 0, 0, 0};
    std::array<std::array<double, Indices::numComponents>, Indices::numPhases> massFraction{};
    massFraction[Indices::Phase::liquid][0] = 0.75;
    massFraction[Indices::Phase::liquid][1] = 0.25;
    massFraction[Indices::Phase::vapor][0] = 0.1;
    massFraction[Indices::Phase::vapor][1] = 0.9;

    const auto result = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Producer,
        2.0,
        90.0,
        pressure,
        density,
        mobility,
        surfaceDensity,
        phaseFraction,
        injectionComposition,
        massFraction);

    near(result.surfacePhaseRate[Indices::Phase::liquid], -2.0, 1e-14,
         "producer oil surface rate");
    near(result.surfacePhaseRate[Indices::Phase::vapor], -200.0, 1e-14,
         "producer gas surface rate");
    near(result.phaseMassRate[Indices::Phase::water], -6000.0, 1e-14,
         "producer water mass rate");
    near(result.componentMassSource[0], -1240.0, 1e-14,
         "component source 0");
    near(result.componentMassSource[1], -760.0, 1e-14,
         "component source 1");

    std::cout << "Well Natural bridge validation: ALL PASS\n";
    return 0;
}
