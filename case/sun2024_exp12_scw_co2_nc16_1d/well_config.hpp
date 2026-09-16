/**
 * @file well_config.hpp
 * @brief Sun-2024 Exp.12 派生算例入口混合注入与出口背压配置。
 */
#pragma once

#include "benchmark_common.hpp"

#include <array>

namespace WellConfig
{
using Type = Sun2024Exp12::Type;
using Control = Sun2024Exp12::Control;
using InjectionPhase = Sun2024Exp12::InjectionPhase;
using Completion = Sun2024Exp12::Completion;

struct WellDefinition : Sun2024Exp12::BaseWellDefinition
{
    // -1 leaves the injector arrays for the case runtime to fill with the
    // common H2O/CO2 mixed-stream definition.
    int injectedComponent{-1};
};

inline static constexpr std::array<WellDefinition, 2> wells{{
    {Sun2024Exp12::baseWells[0], -1},
    {Sun2024Exp12::baseWells[1], -1}
}};
} // namespace WellConfig
