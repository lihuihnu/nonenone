/**
 * @file case_fluid.hpp
 * @brief 超临界水物性与三相相对渗透率闭合。
 */
#pragma once

#include "benchmark_common.hpp"
#include "case_config.hpp"

#include <natural/thermo/aqueous_volume.hpp>

namespace ScwMigrationCaseFluid
{

template <class Indices>
void apply(MPMC::FluidSystem<Indices> &fluid)
{
    using Eval = typename Indices::ValueType;
    ScwMigration3D::applyRelativePermeability(fluid);

    // In legacy independent-water mode B_w is used as rho_res/rho_surface.
    // The IF97 dispatcher covers both the low-density and dense SCW regions.
    fluid.waterFormationVolumeFactor = [](Eval pressure) {
        const Eval reservoirDensity = MPMC::IapwsIf97WaterDensity::density(
            pressure, Eval(ScwMigration3D::temperature));
        return reservoirDensity /
            Eval(CaseConfig::Fluid::surfaceDensity[Indices::Phase::water]);
    };

    // McBride-Wright pure-water limit, evaluated at the fixed case temperature.
    fluid.waterViscosity = [](Eval pressure) {
        constexpr double a = -3.705013;
        constexpr double b = 0.00289258;
        constexpr double c = 3.98950;
        constexpr double d = -0.00326;
        constexpr double t0 = 141.5;
        const Eval reducedPressure = pressure / 1.0e6;
        const double shiftedTemperature =
            ScwMigration3D::temperature / t0 - 1.0;
        return exp(a + b * reducedPressure +
                   (c + d * reducedPressure) / shiftedTemperature) * 1.0e-3;
    };
}

} // namespace ScwMigrationCaseFluid
