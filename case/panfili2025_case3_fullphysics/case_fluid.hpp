#pragma once

#include "case_config.hpp"

#include <common/math.hpp>
#include <natural/fluid_system.hpp>

/**
 * @file case_fluid.hpp
 * @brief Panfili 2025 Case-3 全物理算例的岩石–流体闭合关系。
 *
 * Panfili et al. (2025) show oil-water and gas-oil relative-permeability curves
 * for use case #3 in Fig. 17, but do not publish the underlying numerical table.
 * The Hamilton Eclipse archive should eventually replace these curves.  Until
 * then this file keeps the approximation visible, centralized and testable.
 */
namespace Panfili2025CaseFluid
{

template <class Indices>
void applyRockFluidApproximation(MPMC::FluidSystem<Indices> &fluid)
{
    using Eval = typename Indices::ValueType;

    fluid.waterRelativePermeability = [](Eval sw) {
        constexpr double swc = CaseConfig::LiteratureReference::connateWaterSaturation;
        constexpr double sorw = 0.10; // explicit MPMC closure, not a paper datum
        const double value = MPMC::scalarValue(sw);
        if (value <= swc) return Eval(0.0);
        if (value >= 1.0 - sorw) return Eval(1.0);
        const Eval se = (sw - swc) / (1.0 - swc - sorw);
        return se * se;
    };

    fluid.gasRelativePermeability = [](Eval sg) {
        constexpr double swc = CaseConfig::LiteratureReference::connateWaterSaturation;
        const double value = MPMC::scalarValue(sg);
        if (value <= 0.0) return Eval(0.0);
        if (value >= 1.0 - swc) return Eval(1.0);
        const Eval se = sg / (1.0 - swc);
        return se * se;
    };

    fluid.threePhaseOilRelativePermeability = [](Eval, Eval so, Eval) {
        constexpr double swc = CaseConfig::LiteratureReference::connateWaterSaturation;
        const double value = MPMC::scalarValue(so);
        if (value <= 0.0) return Eval(0.0);
        const Eval se = so / (1.0 - swc);
        return se * se;
    };
}

} // namespace Panfili2025CaseFluid
