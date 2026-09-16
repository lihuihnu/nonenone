/**
 * @file well_config.hpp
 * @brief 传统二维基线算例的注采井位置、控制方式与组分配置。
 */
#pragma once

#include "case_config.hpp"

#include <array>

namespace TraditionalWellConfig
{
using Type = BenchmarkCommon::Type;
using Control = BenchmarkCommon::Control;
using InjectionPhase = BenchmarkCommon::InjectionPhase;
using Completion = BenchmarkCommon::Completion;
struct WellDefinition : BenchmarkCommon::BaseWellDefinition
{
    int injectedComponent{-1};
};
inline static constexpr std::array<WellDefinition, 2> wells{{
    {BenchmarkCommon::baseWells[0], TraditionalConfig::Fluid::co2Component},
    {BenchmarkCommon::baseWells[1], -1}
}};
} // namespace TraditionalWellConfig
