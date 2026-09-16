/**
 * @file case_fluid.hpp
 * @brief 二乘二实验使用的超临界水物性与Corey相对渗透率。
 */
#pragma once

#include "case_config.hpp"

#include <natural/thermo/aqueous_volume.hpp>

namespace Sun2024FactorialFluid
{

template <class Indices>
void apply(MPMC::FluidSystem<Indices> &fluid)
{
    Sun2024Exp12::applyRelativePermeability(fluid);

    using Eval = typename Indices::ValueType;
    fluid.waterFormationVolumeFactor = [](Eval pressure) {
        const Eval reservoirDensity = MPMC::IapwsIf97WaterDensity::density(
            pressure, Eval(Sun2024Factorial::temperature));
        return reservoirDensity /
            Eval(CaseConfig::Fluid::surfaceDensity[Indices::Phase::water]);
    };
    fluid.waterViscosity = [](Eval pressure) {
        constexpr double a = -3.705013;
        constexpr double b = 0.00289258;
        constexpr double c = 3.98950;
        constexpr double d = -0.00326;
        constexpr double t0 = 141.5;
        const Eval reducedPressure = pressure / 1.0e6;
        const double shiftedTemperature =
            Sun2024Factorial::temperature / t0 - 1.0;
        return exp(a + b * reducedPressure +
                   (c + d * reducedPressure) / shiftedTemperature) * 1.0e-3;
    };
}

} // namespace Sun2024FactorialFluid
