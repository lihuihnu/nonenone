#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

/**
 * @file well_config.hpp
 * @brief 五组分优选参数三 EOS 算例的井位、完井和井控参数。
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
    // 西南浅部完井：10 万标准/参考 m3/day 纯 CO2。
    {0, "CO2_INJ", Type::Injector, Control::TotalRate,
     cubicMetresPerDay(100000.0), 66.0 * CaseConfig::bar,
     {1, 4, 0, 2},
     0.10, 0.0, InjectionPhase::Gas, CaseConfig::CommonFluid::co2Component},

    // 东北深部完井：52 bar BHP，位于隔层下方。
    {1, "PROD", Type::Producer, Control::Bhp,
     52.0 * CaseConfig::bar, 52.0 * CaseConfig::bar,
     {CaseConfig::Grid::nx - 2, CaseConfig::Grid::ny - 5, 3, 2},
     0.10, 0.0, InjectionPhase::Oil, -1}
}};

} // namespace WellConfig
