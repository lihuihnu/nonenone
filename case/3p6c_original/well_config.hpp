#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

/**
 * @file well_config.hpp
 * @brief 三相六组分原始对照算例的井位置、完井和控制配置。
 *
 * 新增/删除井只改 wells。i/j=-1 表示该方向最后一个网格；kCount=-1 表示
 * 从 kBegin 完井到顶层。井指数由 main 使用 Peaceman 公式自动计算。
 */
namespace WellConfig
{

using Type = MPMC::cases::well_config::Type;
using Control = MPMC::cases::well_config::Control;
using InjectionPhase = MPMC::cases::well_config::InjectionPhase;
using Completion = MPMC::cases::well_config::Completion<>;
using WellDefinition = MPMC::cases::well_config::StructuredWellDefinition<CaseConfig::Config>;

inline static constexpr std::array<WellDefinition, 2> wells{{
    {0, "INJ", Type::Injector, Control::TotalRate,
     0.008025424544958, 155.0 * CaseConfig::bar,
     {0, 0, 0, -1}, 0.1, 0.0, InjectionPhase::Gas, 1},

    {1, "PROD", Type::Producer, Control::Bhp,
     100.0 * CaseConfig::bar, 100.0 * CaseConfig::bar,
     {-1, -1, 0, -1}, 0.1, 0.0, InjectionPhase::Oil, -1}
}};

} // namespace WellConfig
