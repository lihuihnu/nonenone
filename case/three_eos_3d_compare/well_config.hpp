#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

/**
 * @file well_config.hpp
 * @brief 三 EOS 三维对比算例的井位、完井和井控参数。
 */
namespace WellConfig
{

using Type = MPMC::cases::well_config::Type;
using Control = MPMC::cases::well_config::Control;
using InjectionPhase = MPMC::cases::well_config::InjectionPhase;
using Completion = MPMC::cases::well_config::Completion<>;
using WellDefinition = MPMC::cases::well_config::StructuredWellDefinition<CaseConfig::Config>;

inline constexpr double cubicMetresPerDay(double value)
{
    return value / CaseConfig::secondsPerDay;
}

inline static constexpr std::array<WellDefinition, 2> wells{{
    // 西南下部完井：30 万标准/参考 m3/day 纯 CO2，RATE 控制量按地面体积流量定义。
    {0, "CO2_INJ", Type::Injector, Control::TotalRate,
     cubicMetresPerDay(300000.0), 66.0 * CaseConfig::bar,
     {3, 4, 0, 3},
     0.10, 0.0, InjectionPhase::Gas, CaseConfig::CommonFluid::co2Component},

    // 东北上部完井：50 bar BHP，位于隔夹层上方。
    {1, "PROD", Type::Producer, Control::Bhp,
     50.0 * CaseConfig::bar, 50.0 * CaseConfig::bar,
     {CaseConfig::Grid::nx - 4, CaseConfig::Grid::ny - 5, 4, 4},
     0.10, 0.0, InjectionPhase::Oil, -1}
}};

} // namespace WellConfig
