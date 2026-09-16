/**
 * @file experiment_matrix.hpp
 * @brief 预注册Sun-2024超临界水-CO2因子矩阵与等流量对照。
 */
#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Sun2024Factorial
{

inline constexpr double day = 86400.0;
inline constexpr double minute = 60.0;
inline constexpr double millilitre = 1.0e-6;
inline constexpr double waterCriticalPressure = 22.064e6;
inline constexpr double waterCriticalTemperature = 647.096;
inline constexpr double temperature = 673.15;
inline constexpr double defaultTargetPoreVolumes = 4.0;
inline constexpr double defaultOutputStepPoreVolumes = 0.02;

enum class Role
{
    PhysicalReproduction,
    RateMatchedControl
};

struct ExperimentDefinition
{
    const char *id;
    int sourceExperiment;
    Role role;
    double pressurePa;
    double waterRateMlPerMinute;
    double co2RateMlPerMinute;
    const char *matchedPhysicalExperiment;
};

// The four physical treatments form a 2 x 2 design at 400 degC:
// pressure (23/24 MPa) x added CO2 (0/2 mL/min).  C23/C24 are pure-SCW
// controls whose resolved reservoir rate is matched to R10/R12.  Under the
// in-situ interpretation this is 12 mL/min; under reference-density
// conversion their required reported/reference-equivalent water rate is
// derived rather than assumed.
inline constexpr std::array<ExperimentDefinition, 6> experiments{{
    {"R04", 4, Role::PhysicalReproduction, 23.0e6, 10.0, 0.0, ""},
    {"R06", 6, Role::PhysicalReproduction, 24.0e6, 10.0, 0.0, ""},
    {"R10", 10, Role::PhysicalReproduction, 23.0e6, 10.0, 2.0, ""},
    {"R12", 12, Role::PhysicalReproduction, 24.0e6, 10.0, 2.0, ""},
    {"C23", 0, Role::RateMatchedControl, 23.0e6, 12.0, 0.0, "R10"},
    {"C24", 0, Role::RateMatchedControl, 24.0e6, 12.0, 0.0, "R12"}
}};

enum class RateBasis
{
    InSituVolume,
    ReferenceDensity
};

struct ReferenceDensities
{
    double waterKgPerM3{std::numeric_limits<double>::quiet_NaN()};
    double co2KgPerM3{std::numeric_limits<double>::quiet_NaN()};
};

struct ReservoirDensities
{
    double waterKgPerM3;
    double co2KgPerM3;
};

struct ResolvedInjection
{
    double waterReservoirRateM3PerS;
    double co2ReservoirRateM3PerS;
    double waterMassRateKgPerS;
    double co2MassRateKgPerS;
    double waterReferenceEquivalentMlPerMinute;
    double co2ReferenceEquivalentMlPerMinute;

    [[nodiscard]] double totalReservoirRateM3PerS() const noexcept
    {
        return waterReservoirRateM3PerS + co2ReservoirRateM3PerS;
    }

    [[nodiscard]] double waterReservoirVolumeFraction() const noexcept
    {
        return waterReservoirRateM3PerS / totalReservoirRateM3PerS();
    }

    [[nodiscard]] double co2ReservoirVolumeFraction() const noexcept
    {
        return co2ReservoirRateM3PerS / totalReservoirRateM3PerS();
    }
};

[[nodiscard]] inline const ExperimentDefinition &definition(std::string_view id)
{
    for (const auto &candidate : experiments)
        if (id == candidate.id)
            return candidate;
    throw std::invalid_argument(
        "Unknown -experiment_id '" + std::string(id) +
        "'. Expected R04, R06, R10, R12, C23, or C24.");
}

[[nodiscard]] inline constexpr double toVolumetricRate(double mlPerMinute)
{
    return mlPerMinute * millilitre / minute;
}

[[nodiscard]] inline constexpr double nominalInSituWaterRate(
    const ExperimentDefinition &e)
{
    return toVolumetricRate(e.waterRateMlPerMinute);
}

[[nodiscard]] inline constexpr double nominalInSituCo2Rate(
    const ExperimentDefinition &e)
{
    return toVolumetricRate(e.co2RateMlPerMinute);
}

[[nodiscard]] inline constexpr double nominalInSituTotalRate(
    const ExperimentDefinition &e)
{
    return nominalInSituWaterRate(e) + nominalInSituCo2Rate(e);
}

[[nodiscard]] inline constexpr const char *rateBasisName(RateBasis basis) noexcept
{
    return basis == RateBasis::InSituVolume
        ? "in_situ_volume"
        : "reference_density";
}

[[nodiscard]] inline RateBasis parseRateBasis(std::string_view value)
{
    if (value == "in_situ_volume")
        return RateBasis::InSituVolume;
    if (value == "reference_density")
        return RateBasis::ReferenceDensity;
    throw std::invalid_argument(
        "Unknown -rate_basis '" + std::string(value) +
        "'. Expected in_situ_volume or reference_density.");
}

inline void validateDensities(
    const ReferenceDensities &reference,
    const ReservoirDensities &reservoir,
    bool requireCo2Reference)
{
    if (!(reservoir.waterKgPerM3 > 0.0) ||
        !(reservoir.co2KgPerM3 > 0.0) ||
        !std::isfinite(reservoir.waterKgPerM3) ||
        !std::isfinite(reservoir.co2KgPerM3))
        throw std::invalid_argument(
            "Reservoir phase densities must be finite and positive.");
    if (!(reference.waterKgPerM3 > 0.0) ||
        !std::isfinite(reference.waterKgPerM3))
        throw std::invalid_argument(
            "reference_density requires -water_reference_density_kg_m3.");
    if (requireCo2Reference &&
        (!(reference.co2KgPerM3 > 0.0) ||
         !std::isfinite(reference.co2KgPerM3)))
        throw std::invalid_argument(
            "This run requires -co2_reference_density_kg_m3.");
}

[[nodiscard]] inline ResolvedInjection resolvePhysicalInjection(
    const ExperimentDefinition &e,
    RateBasis basis,
    const ReferenceDensities &reference,
    const ReservoirDensities &reservoir)
{
    const double reportedWater = toVolumetricRate(e.waterRateMlPerMinute);
    const double reportedCo2 = toVolumetricRate(e.co2RateMlPerMinute);
    if (basis == RateBasis::ReferenceDensity)
        validateDensities(reference, reservoir, reportedCo2 > 0.0);

    const double waterMass = basis == RateBasis::InSituVolume
        ? reportedWater * reservoir.waterKgPerM3
        : reportedWater * reference.waterKgPerM3;
    const double co2Mass = basis == RateBasis::InSituVolume
        ? reportedCo2 * reservoir.co2KgPerM3
        : reportedCo2 * reference.co2KgPerM3;

    return {
        waterMass / reservoir.waterKgPerM3,
        co2Mass / reservoir.co2KgPerM3,
        waterMass,
        co2Mass,
        e.waterRateMlPerMinute,
        e.co2RateMlPerMinute};
}

[[nodiscard]] inline ResolvedInjection resolveInjection(
    const ExperimentDefinition &e,
    RateBasis basis,
    const ReferenceDensities &reference,
    const ReservoirDensities &reservoir)
{
    if (e.role == Role::PhysicalReproduction)
        return resolvePhysicalInjection(e, basis, reference, reservoir);

    const auto &matched = definition(e.matchedPhysicalExperiment);
    const auto physical = resolvePhysicalInjection(
        matched, basis, reference, reservoir);
    const double waterReservoirRate = physical.totalReservoirRateM3PerS();
    const double waterMass = waterReservoirRate * reservoir.waterKgPerM3;
    const double equivalentReferenceRate = basis == RateBasis::InSituVolume
        ? waterReservoirRate
        : waterMass / reference.waterKgPerM3;
    return {
        waterReservoirRate,
        0.0,
        waterMass,
        0.0,
        equivalentReferenceRate / millilitre * minute,
        0.0};
}

[[nodiscard]] inline double durationDays(
    const ResolvedInjection &injection,
    double poreVolumeM3,
    double injectedPoreVolumes)
{
    if (!(poreVolumeM3 > 0.0) || !(injectedPoreVolumes > 0.0))
        throw std::invalid_argument("Pore volume and target PVI must be positive.");
    return injectedPoreVolumes * poreVolumeM3 /
        injection.totalReservoirRateM3PerS() / day;
}

[[nodiscard]] inline const char *roleName(Role role) noexcept
{
    return role == Role::PhysicalReproduction
        ? "physical_reproduction"
        : "rate_matched_control";
}

} // namespace Sun2024Factorial
