#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>

/**
 * @file well_config.hpp
 * @brief Panfili 2025 Case-3 全物理算例的井位置、完井和控制配置。
 *
 * The paper states that the SAME two crest wells produce for 32 years, remain
 * idle for two years, and are then converted to CO2 injectors.  v25 represented
 * the producer and injector stages as four separate Well objects sharing two
 * BHP representative cells.  Natural intentionally forbids that layout because
 * each representative cell owns exactly one well-pressure equation.
 *
 * v26 therefore keeps exactly two physical wells and changes their operating
 * mode with time through NaturalPetscRuntime::setWellScheduleUpdater().  This
 * matches the paper semantics and preserves one BHP unknown per physical well.
 * Exact Hamilton completion-cell indices remain unavailable in the paper; the
 * structured average-property wrapper keeps the same two symmetric crest cells.
 */
namespace WellConfig
{

using Type = MPMC::cases::well_config::Type;
using Control = MPMC::cases::well_config::Control;
using InjectionPhase = MPMC::cases::well_config::InjectionPhase;
using Completion = MPMC::cases::well_config::Completion<1>;

enum class OperatingStage { Depletion, Idle, Injection, Closed };

struct WellDefaults
{
    static constexpr Control control = Control::GasRate;
    inline static constexpr Completion completion{0, 0, 0, 1};
    static constexpr int injectedComponent = CaseConfig::Fluid::co2Component;
    static constexpr double minimumBhp = CaseConfig::LiteratureReference::producerMinimumBhp;
    static constexpr double maximumBhp = CaseConfig::LiteratureReference::injectorMaximumBhp;
    static constexpr bool scheduleEnabled = true;
    static constexpr double openTime = 0.0;
    static constexpr double closeTime =
        CaseConfig::LiteratureReference::figure25InjectionEndYear *
        CaseConfig::daysPerYear * CaseConfig::secondsPerDay;
};

using WellDefinition =
    MPMC::cases::well_config::ScheduledStructuredWellDefinition<
        CaseConfig::Config, Completion, WellDefaults>;

inline constexpr double years(double value)
{
    return value * CaseConfig::daysPerYear * CaseConfig::secondsPerDay;
}

inline constexpr double mscfPerDayToM3PerSecond(double value)
{
    // Petroleum Mscf = 1000 standard cubic feet.
    return value * 1000.0 * CaseConfig::standardCubicFoot /
           CaseConfig::secondsPerDay;
}


inline constexpr double producerTarget()
{
    return mscfPerDayToM3PerSecond(
        CaseConfig::LiteratureReference::producerGasRateMscfPerDay);
}

inline constexpr double injectorTarget()
{
    return mscfPerDayToM3PerSecond(
        CaseConfig::LiteratureReference::injectorGasRateMscfPerDay);
}

inline constexpr OperatingStage stageAt(double time)
{
    // Backward-Euler attempts are evaluated at the attempted end time.  Keeping
    // t=32 y in depletion makes the 31->32 y step a production step; keeping
    // t=34 y idle makes the first injection step 34->35 y.
    if (time <= years(CaseConfig::LiteratureReference::depletionYears))
        return OperatingStage::Depletion;
    if (time <= years(CaseConfig::LiteratureReference::injectionStartYear))
        return OperatingStage::Idle;
    if (time <= years(CaseConfig::LiteratureReference::figure25InjectionEndYear))
        return OperatingStage::Injection;
    return OperatingStage::Closed;
}

inline constexpr Completion crestA{
    CaseConfig::Grid::nx / 3,
    CaseConfig::Grid::ny / 2,
    CaseConfig::Grid::nz - 1,
    1};
inline constexpr Completion crestB{
    2 * CaseConfig::Grid::nx / 3,
    CaseConfig::Grid::ny / 2,
    CaseConfig::Grid::nz - 1,
    1};

inline static constexpr std::array<WellDefinition, 2> wells{{
    {0, "CREST_A", Type::Producer, Control::GasRate, producerTarget(),
     CaseConfig::InitialState::pressure, crestA, 0.10, 0.0,
     InjectionPhase::Gas, CaseConfig::Fluid::co2Component,
     CaseConfig::LiteratureReference::producerMinimumBhp,
     CaseConfig::LiteratureReference::injectorMaximumBhp,
     true, 0.0, years(CaseConfig::LiteratureReference::figure25InjectionEndYear)},
    {1, "CREST_B", Type::Producer, Control::GasRate, producerTarget(),
     CaseConfig::InitialState::pressure, crestB, 0.10, 0.0,
     InjectionPhase::Gas, CaseConfig::Fluid::co2Component,
     CaseConfig::LiteratureReference::producerMinimumBhp,
     CaseConfig::LiteratureReference::injectorMaximumBhp,
     true, 0.0, years(CaseConfig::LiteratureReference::figure25InjectionEndYear)}
}};

} // namespace WellConfig
