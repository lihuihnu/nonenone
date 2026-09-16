#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

/**
 * @file well_config.hpp
 * @brief Ma 2021 五组分油气水三相文献算例的井位置、完井和控制配置。
 *
 * Ma et al. (2021) provide a thermodynamic flash benchmark, not reservoir wells.
 * These two deliberately mild wells are therefore MPMC-side transport test
 * conditions.  A pure-CH4 gas injector exercises hydrocarbon transfer into the
 * water-rich phase; a small BHP drawdown produces the five-component mixture.
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
    {0, "C1_INJ", Type::Injector, Control::TotalRate,
     cubicMetresPerDay(1.0), 14.0 * CaseConfig::bar,
     {0, 0, 0, 1}, 0.10, 0.0,
     InjectionPhase::Gas, CaseConfig::Fluid::methaneComponent},

    {1, "PROD", Type::Producer, Control::Bhp,
     13.50 * CaseConfig::bar, 13.50 * CaseConfig::bar,
     {-1, 0, 0, 1}, 0.10, 0.0,
     InjectionPhase::Oil, -1}
}};

} // namespace WellConfig
