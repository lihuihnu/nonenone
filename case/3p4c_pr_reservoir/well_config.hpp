#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

/**
 * @file well_config.hpp
 * @brief H2O–CO2–CH4–nC16 三相四组分 PR 储层算例的井位置、完井和控制配置。
 *
 * The injector supplies pure CO2 through the gas injection slot.  The producer
 * uses a mild BHP drawdown so the first integration test exercises transport
 * without intentionally forcing a violent phase-boundary crossing in one step.
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
    // 100 m3/day pure CO2 at the left-center cell.
    {0, "CO2_INJ", Type::Injector, Control::TotalRate,
     cubicMetresPerDay(100.0), 205.0 * CaseConfig::bar,
     {0, CaseConfig::Grid::ny / 2, 0, 1},
     0.10, 0.0, InjectionPhase::Gas, CaseConfig::Fluid::co2Component},

    // Mild 5 bar drawdown at the right-center cell.
    {1, "PROD", Type::Producer, Control::Bhp,
     195.0 * CaseConfig::bar, 195.0 * CaseConfig::bar,
     {-1, CaseConfig::Grid::ny / 2, 0, 1},
     0.10, 0.0, InjectionPhase::Oil, -1}
}};

} // namespace WellConfig
