/**
 * @file sun2024_factorial_case_test.cpp
 * @brief 可复现Sun-2024实验矩阵的无PETSc契约测试。
 */
#include "../../../case/sun2024_scw_co2_nc16_factorial_1d/case_config.hpp"
#include "../../../case/sun2024_scw_co2_nc16_factorial_1d/case_fluid.hpp"
#include "../../../case/sun2024_scw_co2_nc16_factorial_1d/experiment_matrix.hpp"

#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/aqueous_volume.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance,
          const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(message);
}

void checkMatrix()
{
    static_assert(!Indices::fullyCompositionalThreePhase);
    static_assert(Indices::hasIndependentWaterConservation);
    static_assert(Indices::numComponents == 2);
    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();

    require(Sun2024Factorial::experiments.size() == 6,
            "factorial must contain four physical runs and two controls");
    const auto &r04 = Sun2024Factorial::definition("R04");
    const auto &r06 = Sun2024Factorial::definition("R06");
    const auto &r10 = Sun2024Factorial::definition("R10");
    const auto &r12 = Sun2024Factorial::definition("R12");
    const auto &c23 = Sun2024Factorial::definition("C23");
    const auto &c24 = Sun2024Factorial::definition("C24");

    near(r04.pressurePa, r10.pressurePa, 0.0, "23 MPa block");
    near(r06.pressurePa, r12.pressurePa, 0.0, "24 MPa block");
    near(r04.waterRateMlPerMinute, 10.0, 0.0, "R04 water rate");
    near(r10.co2RateMlPerMinute, 2.0, 0.0, "R10 CO2 rate");
    const Sun2024Factorial::ReservoirDensities reservoir{100.0, 200.0};
    const Sun2024Factorial::ReferenceDensities reference{1000.0, 800.0};
    const auto r12InSitu = Sun2024Factorial::resolveInjection(
        r12, Sun2024Factorial::RateBasis::InSituVolume, {}, reservoir);
    const auto c24InSitu = Sun2024Factorial::resolveInjection(
        c24, Sun2024Factorial::RateBasis::InSituVolume, {}, reservoir);
    const auto r10InSitu = Sun2024Factorial::resolveInjection(
        r10, Sun2024Factorial::RateBasis::InSituVolume, {}, reservoir);
    const auto c23InSitu = Sun2024Factorial::resolveInjection(
        c23, Sun2024Factorial::RateBasis::InSituVolume, {}, reservoir);
    near(r10InSitu.totalReservoirRateM3PerS(),
         c23InSitu.totalReservoirRateM3PerS(), 1.0e-14,
         "C23 must match R10 in-situ total rate");
    near(r12InSitu.totalReservoirRateM3PerS(),
         c24InSitu.totalReservoirRateM3PerS(), 1.0e-14,
         "C24 must match R12 in-situ total rate");
    near(r12InSitu.waterReservoirVolumeFraction(), 5.0 / 6.0,
         1.0e-14, "R12 in-situ SCW fraction");
    near(r12InSitu.co2ReservoirVolumeFraction(), 1.0 / 6.0,
         1.0e-14, "R12 in-situ CO2 fraction");
    near(c24InSitu.waterReservoirVolumeFraction(), 1.0,
         1.0e-14, "C24 pure SCW fraction");

    const auto r12Reference = Sun2024Factorial::resolveInjection(
        r12, Sun2024Factorial::RateBasis::ReferenceDensity,
        reference, reservoir);
    const auto c24Reference = Sun2024Factorial::resolveInjection(
        c24, Sun2024Factorial::RateBasis::ReferenceDensity,
        reference, reservoir);
    near(r12Reference.waterReservoirRateM3PerS,
         Sun2024Factorial::toVolumetricRate(100.0), 1.0e-14,
         "water reference-density conversion");
    near(r12Reference.co2ReservoirRateM3PerS,
         Sun2024Factorial::toVolumetricRate(8.0), 1.0e-14,
         "CO2 reference-density conversion");
    near(r12Reference.totalReservoirRateM3PerS(),
         c24Reference.totalReservoirRateM3PerS(), 1.0e-14,
         "C24 must match converted R12 reservoir rate");
    near(c24Reference.waterReferenceEquivalentMlPerMinute, 10.8,
         1.0e-14, "C24 derived reference-state water rate");

    const double expectedPv = 0.48 *
        (3.14159265358979323846 * 0.039 * 0.039 / 4.0) * 0.39;
    near(Sun2024Exp12::poreVolume, expectedPv, 1.0e-14,
         "sandpack pore volume");
    near(Sun2024Factorial::durationDays(r12InSitu, expectedPv, 1.0) *
             Sun2024Factorial::day / 60.0,
         18.6356134618293, 1.0e-12, "R12 minutes per PVI");
    require(r04.pressurePa > Sun2024Factorial::waterCriticalPressure,
            "all registered pressures must exceed the water critical point");
}

void checkFluid()
{
    auto fluid = MPMC::cases::makeFluidSystem<
        Indices, CaseConfig::PrFactoryConfig>();
    Sun2024FactorialFluid::apply(fluid);

    for (const auto &experiment : Sun2024Factorial::experiments)
    {
        const double expectedDensity = MPMC::IapwsIf97WaterDensity::density(
            experiment.pressurePa, Sun2024Factorial::temperature);
        const double configuredDensity =
            fluid.waterFormationVolumeFactor(experiment.pressurePa) *
            CaseConfig::Fluid::surfaceDensity[Indices::Phase::water];
        near(configuredDensity, expectedDensity, 1.0e-12,
             std::string(experiment.id) + " IF97 density");
        const double viscosity = fluid.waterViscosity(experiment.pressurePa);
        require(std::isfinite(viscosity) && viscosity > 0.0 && viscosity < 2.0e-4,
                std::string(experiment.id) + " SCW viscosity");
    }
}

} // namespace

int main()
{
    try
    {
        checkMatrix();
        checkFluid();
        std::cout << "[PASS] sun2024_factorial_case_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] sun2024_factorial_case_test: "
                  << error.what() << '\n';
        return 1;
    }
}
