/**
 * @file well_config.hpp
 * @brief 基础井对象定义，实验流量与压力目标在运行时赋值。
 */
#pragma once

#include "../sun2024_exp12_scw_co2_nc16_1d/benchmark_common.hpp"

#include <array>

namespace WellConfig
{
using Type = Sun2024Exp12::Type;
using Control = Sun2024Exp12::Control;
using InjectionPhase = Sun2024Exp12::InjectionPhase;
using Completion = Sun2024Exp12::Completion;

struct WellDefinition : Sun2024Exp12::BaseWellDefinition
{
    int injectedComponent{-1};
};

inline static constexpr std::array<WellDefinition, 2> wells{{
    {Sun2024Exp12::baseWells[0], -1},
    {Sun2024Exp12::baseWells[1], -1}
}};
} // namespace WellConfig
