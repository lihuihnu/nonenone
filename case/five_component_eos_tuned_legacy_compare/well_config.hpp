/**
 * @file well_config.hpp
 * @brief 五组分传统独立水 PR 对照算例的井定义与控制参数。
 */
#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

namespace LegacyWellConfig
{

using Type = MPMC::cases::well_config::Type;
using Control = MPMC::cases::well_config::Control;
using InjectionPhase = MPMC::cases::well_config::InjectionPhase;
using Completion = MPMC::cases::well_config::Completion<>;
using WellDefinition = MPMC::cases::well_config::StructuredWellDefinition<LegacyCaseConfig::Config>;

inline constexpr double cubicMetresPerDay(double value)
{
    return value / LegacyCaseConfig::secondsPerDay;
}

inline static constexpr std::array<WellDefinition, 2> wells{{
    {0, "CO2_INJ", Type::Injector, Control::TotalRate,
     cubicMetresPerDay(100000.0), 66.0 * LegacyCaseConfig::bar,
     {1, 4, 0, 2}, 0.10, 0.0, InjectionPhase::Gas,
     LegacyCaseConfig::CommonFluid::co2Component},
    {1, "PROD", Type::Producer, Control::Bhp,
     52.0 * LegacyCaseConfig::bar, 52.0 * LegacyCaseConfig::bar,
     {LegacyCaseConfig::Grid::nx - 2, LegacyCaseConfig::Grid::ny - 5, 3, 2},
     0.10, 0.0, InjectionPhase::Oil, -1}
}};

} // namespace LegacyWellConfig
