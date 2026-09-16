/**
 * @file case_validation.hpp
 * @brief 不依赖 PETSc 的算例配置合法性检查。
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace MPMC::cases
{

/**
 * @brief 在不依赖 PETSc 的条件下校验用户侧算例配置。
 *
 * Fully-compositional O/G/W cases intentionally expose a different initialization
 * contract from legacy cases.  They provide P-T-z only and must not be forced to
 * invent phase saturations/compositions in `case_config.hpp`; those quantities
 * are flash results.  Keeping this validation in a PETSc-free header lets unit
 * tests validate the exact production case configuration.
 */
template <class Indices, class Config>
void validateCaseConfig()
{
    const auto finite = [](double v) { return std::isfinite(v); };

    if (!(Config::InitialState::pressure > 0.0) || !finite(Config::InitialState::pressure))
        throw std::invalid_argument("Initial pressure must be finite and positive.");

    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        if (!(Config::InitialState::temperature > 0.0) ||
            !finite(Config::InitialState::temperature))
            throw std::invalid_argument("Initial temperature must be finite and positive.");

        const double temperatureTolerance =
            1.0e-10 * std::max(1.0, std::abs(Config::Fluid::temperature));
        if (std::abs(Config::InitialState::temperature - Config::Fluid::temperature) >
            temperatureTolerance)
            throw std::invalid_argument(
                "The current Natural model is isothermal: InitialState::temperature must match Fluid::temperature.");

        if (Config::Fluid::waterComponent < 0 ||
            Config::Fluid::waterComponent >= Indices::numComponents)
            throw std::invalid_argument("Fluid::waterComponent is outside the component range.");

        double total = 0.0;
        for (double value : Config::InitialState::overallComposition)
        {
            if (value < 0.0 || value > 1.0 || !finite(value))
                throw std::invalid_argument(
                    "Initial overall composition contains an invalid mole fraction.");
            total += value;
        }
        if (std::abs(total - 1.0) > 1.0e-8)
            throw std::invalid_argument(
                "Initial overall composition must sum to one in the full three-phase mode.");
    }
    else
    {
        const std::array<double, 3> saturationValues{
            Config::InitialState::oilSaturation,
            Config::InitialState::gasSaturation,
            Indices::hasWater ? Config::InitialState::waterSaturation : 0.0};
        for (double value : saturationValues)
            if (value < 0.0 || value > 1.0 || !finite(value))
                throw std::invalid_argument("Initial saturation must be finite and in [0,1].");

        if constexpr (!Config::InitialState::preserveReferenceValues)
        {
            double saturation =
                Config::InitialState::oilSaturation + Config::InitialState::gasSaturation;
            if constexpr (Indices::hasWater)
                saturation += Config::InitialState::waterSaturation;
            if (std::abs(saturation - 1.0) > 1.0e-6)
                throw std::invalid_argument("Initial phase saturations must sum to 1.");
        }

        const auto checkComposition = [](const auto &composition, const char *label) {
            double sum = 0.0;
            for (double value : composition)
            {
                if (value < 0.0 || value > 1.0 || !std::isfinite(value))
                    throw std::invalid_argument(std::string(label) + " contains an invalid value.");
                sum += value;
            }
            if constexpr (!Config::InitialState::preserveReferenceValues)
                if (std::abs(sum - 1.0) > 1.0e-6)
                    throw std::invalid_argument(std::string(label) + " must sum to 1.");
        };
        checkComposition(Config::InitialState::oilComposition, "Initial oil composition");
        checkComposition(Config::InitialState::gasComposition, "Initial gas composition");
    }

    if constexpr (Indices::hasAqueousCO2Dissolution)
    {
        if (Config::Dissolution::component < 0 ||
            Config::Dissolution::component >= Indices::numComponents)
            throw std::invalid_argument("Dissolution component is outside the component range.");
        if (!(Config::Dissolution::waterMolarMass > 0.0))
            throw std::invalid_argument("Dissolution water molar mass must be positive.");
        if (Config::Dissolution::initialWaterCO2MoleFraction < 0.0 ||
            Config::Dissolution::initialWaterCO2MoleFraction >= 1.0)
            throw std::invalid_argument("Initial aqueous CO2 mole fraction must be in [0,1). ");
    }

    if constexpr (Indices::hasLandTrapping)
    {
        if (!(Config::Land::constant > 0.0) || !finite(Config::Land::constant))
            throw std::invalid_argument("Land constant must be finite and positive.");
    }

    if constexpr (Indices::hasAdsorption)
    {
        if (!(Config::Adsorption::rockDensity > 0.0))
            throw std::invalid_argument("Adsorption rock density must be positive.");
        for (int c = 0; c < Indices::numComponents; ++c)
        {
            const auto i = static_cast<std::size_t>(c);
            if (Config::Adsorption::thetaMax[i] < 0.0 ||
                Config::Adsorption::coefficient[i] < 0.0)
                throw std::invalid_argument("Adsorption parameters must be non-negative.");
        }
    }
}

} // namespace MPMC::cases
