/**
 * @file scw_co2_nc4_nc16_3d_case_test.cpp
 * @brief 三维 SCW/CO2/nC4/nC16 算例的 PETSc-free 配置与物性预检。
 */
#include "../../../case/scw_co2_nc4_nc16_3d_migration/case_config.hpp"
#include "../../../case/scw_co2_nc4_nc16_3d_migration/case_fluid.hpp"
#include "../../../case/scw_co2_nc4_nc16_3d_migration/well_config.hpp"

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

void checkConfiguration()
{
    static_assert(!Indices::fullyCompositionalThreePhase);
    static_assert(Indices::hasIndependentWaterConservation);
    static_assert(Indices::numComponents == 3);
    static_assert(CaseConfig::Grid::nx * CaseConfig::Grid::ny *
                      CaseConfig::Grid::nz == 2304);

    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();
    require(ScwMigration3D::temperature >
                ScwMigration3D::waterCriticalTemperature,
            "temperature must exceed the water critical point");
    require(ScwMigration3D::initialPressure >
                ScwMigration3D::waterCriticalPressure,
            "pressure must exceed the water critical point");
    near(CaseConfig::InitialState::oilSaturation +
             CaseConfig::InitialState::gasSaturation +
             CaseConfig::InitialState::waterSaturation,
         1.0, 1.0e-13, "initial saturations");
    near(CaseConfig::InitialState::oilComposition[0] +
             CaseConfig::InitialState::oilComposition[1] +
             CaseConfig::InitialState::oilComposition[2],
         1.0, 1.0e-13, "initial oil composition");
    near(ScwMigration3D::totalInjectionRate * ScwMigration3D::year /
             ScwMigration3D::nominalPoreVolume,
         ScwMigration3D::injectionPoreVolumesPerYear, 1.0e-13,
         "injection PV rate");
    near(ScwMigration3D::defaultScwInjectionVolumeFraction,
         0.75, 1.0e-13, "default SCW fraction");

    require(WellConfig::wells.size() == 2,
            "case must have one injector and one producer");
    near(WellConfig::wells[0].target,
         ScwMigration3D::totalInjectionRate, 1.0e-13,
         "injector reservoir total rate");
    near(WellConfig::wells[1].target,
         ScwMigration3D::totalInjectionRate, 1.0e-13,
         "producer must balance the reservoir total injection rate");
    require(WellConfig::wells[0].completion.kBegin <
                WellConfig::wells[1].completion.kBegin,
            "injector must be completed below producer");
    require(CaseConfig::Rock::kz(12, CaseConfig::Rock::channelCenter(12), 2) >
                CaseConfig::Rock::kz(5, CaseConfig::Rock::channelCenter(5), 2),
            "baffle must retain a localized vertical-flow window");
    require(CaseConfig::Rock::kx(2, CaseConfig::Rock::channelCenter(2), 0) >
                CaseConfig::Rock::kx(2, 15, 0),
            "injector must lie in a high-permeability channel");
}

void checkProductionFluidClosures()
{
    auto fluid = MPMC::cases::makeFluidSystem<
        Indices, CaseConfig::PrFactoryConfig>();
    ScwMigrationCaseFluid::apply(fluid);

    static_assert(CaseConfig::Fluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::PengRobinson,
                  "hydrocarbon subsystem must use Peng-Robinson");
    const double density = MPMC::IapwsIf97WaterDensity::density(
        ScwMigration3D::initialPressure, ScwMigration3D::temperature);
    const double recoveredDensity =
        fluid.waterFormationVolumeFactor(ScwMigration3D::initialPressure) *
        CaseConfig::Fluid::surfaceDensity[Indices::Phase::water];
    require(std::isfinite(density) && density > 0.0,
            "SCW density must be finite and positive");
    near(recoveredDensity, density, 1.0e-12,
         "case must use the production IF97 density closure");

    const double viscosity =
        fluid.waterViscosity(ScwMigration3D::initialPressure);
    require(std::isfinite(viscosity) && viscosity > 0.0 && viscosity < 2.0e-4,
            "SCW viscosity must be finite, positive and in the expected range");
    require(fluid.waterRelativePermeability(0.30) > 0.0,
            "mobile SCW must have positive relative permeability");
    require(fluid.gasRelativePermeability(0.20) > 0.0,
            "CO2-rich gas must have positive relative permeability");
    require(fluid.threePhaseOilRelativePermeability(0.15, 0.70, 0.15) > 0.0,
            "oil must remain mobile away from residual saturation");
}

} // namespace

int main()
{
    try
    {
        checkConfiguration();
        checkProductionFluidClosures();
        std::cout << "[PASS] scw_co2_nc4_nc16_3d_case_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] scw_co2_nc4_nc16_3d_case_test: "
                  << error.what() << '\n';
        return 1;
    }
}
