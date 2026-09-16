/**
 * @file indices.hpp
 * @brief Natural 方程、主变量和相状态的编译期索引布局。
 */
#pragma once

#include <indices/model_config.hpp>

#include <ad/Evaluation.hpp>

#include <array>
#include <cstddef>
#include <type_traits>

namespace MPMC
{
namespace detail
{
template <std::size_t Size>
[[nodiscard]] constexpr std::array<int, Size>
makeIndexArray(int offset) noexcept
{
    std::array<int, Size> result{};
    for (std::size_t i = 0; i < Size; ++i)
        result[i] = offset + static_cast<int>(i);
    return result;
}
} // namespace detail

/**
 * @brief Natural formulation 的编译期变量与方程索引布局。
 *
 * Legacy mode deliberately preserves the v16 numeric layout.  The new fully
 * compositional three-phase mode uses
 *
 *   p,
 *   x_o[0..N-2], x_g[0..N-2], x_w[0..N-2],
 *   S_o, S_g, S_w,
 *   BHP(optional)
 *
 * as primary variables.  Its equations are N component balances, N oil-gas
 * fugacity equalities, N oil-water fugacity equalities, one saturation closure
 * and an optional well equation.  Thus the local nonlinear system remains
 * square for arbitrary N.
 */
template <class Config, bool UseAutomaticDifferentiation>
struct Indices final
{
    using ModelConfig = Config;
    static constexpr int disabledIndex = -1;

    static constexpr int numComponents = Config::numComponents;
    static constexpr bool hasWater = Config::hasWater;
    static constexpr bool hasWellUnknown = Config::hasWellUnknown;
    static constexpr bool hasAqueousCO2Dissolution = Config::hasAqueousCO2Dissolution;
    static constexpr bool hasAdsorption = Config::hasAdsorption;
    static constexpr bool hasLandTrapping = Config::hasLandTrapping;
    static constexpr bool fullyCompositionalThreePhase = Config::fullyCompositionalThreePhase;
    static constexpr bool hasIndependentWaterConservation = Config::hasIndependentWaterConservation;
    static constexpr int numHydrocarbonPhases = Config::numHydrocarbonPhases;
    static constexpr int numThermodynamicPhases = Config::numThermodynamicPhases;
    static constexpr int numPhases = Config::numPhases;
    static constexpr bool useAutomaticDifferentiation = UseAutomaticDifferentiation;

    static constexpr int numIndependentCompositionsPerPhase = numComponents - 1;

    struct Phase final
    {
        static constexpr int liquid = 0;   ///< oil-rich liquid
        static constexpr int vapor = 1;    ///< vapor/gas
        static constexpr int water = hasWater ? 2 : disabledIndex; ///< water-rich liquid in full mode
        static constexpr int count = numPhases;
    };

    struct Primary final
    {
        static constexpr int pressure = 0;
        static constexpr int liquidCompositionBegin = pressure + 1;
        static constexpr auto liquidComposition =
            detail::makeIndexArray<static_cast<std::size_t>(numIndependentCompositionsPerPhase)>(
                liquidCompositionBegin);
        static constexpr int afterLiquidComposition =
            liquidCompositionBegin + numIndependentCompositionsPerPhase;

        // ------------------------ fully compositional ---------------------
        static constexpr int fullVaporCompositionBegin = afterLiquidComposition;
        static constexpr int fullAfterVaporComposition =
            fullVaporCompositionBegin + numIndependentCompositionsPerPhase;
        static constexpr int aqueousCompositionBegin =
            fullyCompositionalThreePhase ? fullAfterVaporComposition : disabledIndex;
        static constexpr auto aqueousComposition =
            detail::makeIndexArray<static_cast<std::size_t>(numIndependentCompositionsPerPhase)>(
                fullyCompositionalThreePhase ? aqueousCompositionBegin : 0);
        // Preferred full-model naming.  `aqueousComposition` remains as a
        // compatibility alias for early v17 development code.
        static constexpr int waterCompositionBegin = aqueousCompositionBegin;
        static constexpr auto waterComposition = aqueousComposition;
        static constexpr int fullAfterAqueousComposition =
            fullyCompositionalThreePhase
                ? waterCompositionBegin + numIndependentCompositionsPerPhase
                : disabledIndex;

        // ----------------------------- legacy ----------------------------
        static constexpr int legacyWaterSaturation =
            hasWater && !fullyCompositionalThreePhase
                ? afterLiquidComposition
                : disabledIndex;
        static constexpr int legacyAfterWaterSaturation =
            afterLiquidComposition + ((hasWater && !fullyCompositionalThreePhase) ? 1 : 0);
        static constexpr int legacyWellPressure =
            !fullyCompositionalThreePhase && hasWellUnknown
                ? legacyAfterWaterSaturation
                : disabledIndex;
        static constexpr int legacyAfterWellPressure =
            legacyAfterWaterSaturation + ((!fullyCompositionalThreePhase && hasWellUnknown) ? 1 : 0);
        static constexpr int legacyLiquidSaturation = legacyAfterWellPressure;
        static constexpr int legacyVaporCompositionBegin = legacyLiquidSaturation + 1;
        static constexpr int legacyVaporSaturation =
            legacyVaporCompositionBegin + numIndependentCompositionsPerPhase;
        static constexpr int legacyBaseCount = legacyVaporSaturation + 1;
        static constexpr int legacyAqueousCO2 =
            hasAqueousCO2Dissolution ? legacyBaseCount : disabledIndex;
        static constexpr int legacyCount =
            legacyBaseCount + (hasAqueousCO2Dissolution ? 1 : 0);

        // --------------------- selected public indices -------------------
        static constexpr int vaporCompositionBegin =
            fullyCompositionalThreePhase
                ? fullVaporCompositionBegin
                : legacyVaporCompositionBegin;
        static constexpr auto vaporComposition =
            detail::makeIndexArray<static_cast<std::size_t>(numIndependentCompositionsPerPhase)>(
                vaporCompositionBegin);

        static constexpr int liquidSaturation =
            fullyCompositionalThreePhase
                ? fullAfterAqueousComposition
                : legacyLiquidSaturation;
        static constexpr int vaporSaturation =
            fullyCompositionalThreePhase
                ? liquidSaturation + 1
                : legacyVaporSaturation;
        static constexpr int waterSaturation =
            fullyCompositionalThreePhase
                ? vaporSaturation + 1
                : legacyWaterSaturation;

        static constexpr int fullAfterSaturations =
            fullyCompositionalThreePhase ? waterSaturation + 1 : disabledIndex;
        static constexpr int wellPressure =
            fullyCompositionalThreePhase
                ? (hasWellUnknown ? fullAfterSaturations : disabledIndex)
                : legacyWellPressure;
        static constexpr int fullCount =
            fullyCompositionalThreePhase
                ? fullAfterSaturations + (hasWellUnknown ? 1 : 0)
                : disabledIndex;

        static constexpr int baseCount =
            fullyCompositionalThreePhase ? fullCount : legacyBaseCount;
        static constexpr int aqueousCO2MoleFraction =
            fullyCompositionalThreePhase ? disabledIndex : legacyAqueousCO2;
        static constexpr int count =
            fullyCompositionalThreePhase ? fullCount : legacyCount;
    };

    struct Equation final
    {
        static constexpr int massConservationBegin = 0;
        static constexpr auto massConservation =
            detail::makeIndexArray<static_cast<std::size_t>(numComponents)>(massConservationBegin);
        static constexpr int afterMassConservation = numComponents;

        // Legacy water / well placement is preserved numerically.
        static constexpr int waterConservation =
            hasIndependentWaterConservation ? afterMassConservation : disabledIndex;
        static constexpr int legacyAfterWaterConservation =
            afterMassConservation + (hasIndependentWaterConservation ? 1 : 0);
        static constexpr int legacyWellControl =
            !fullyCompositionalThreePhase && hasWellUnknown
                ? legacyAfterWaterConservation
                : disabledIndex;
        static constexpr int legacyAfterWellControl =
            legacyAfterWaterConservation + ((!fullyCompositionalThreePhase && hasWellUnknown) ? 1 : 0);

        static constexpr int fugacityBegin =
            fullyCompositionalThreePhase
                ? afterMassConservation
                : legacyAfterWellControl;
        static constexpr auto fugacity =
            detail::makeIndexArray<static_cast<std::size_t>(numComponents)>(fugacityBegin);

        static constexpr int aqueousFugacityBegin =
            fullyCompositionalThreePhase ? fugacityBegin + numComponents : disabledIndex;
        static constexpr auto aqueousFugacity =
            detail::makeIndexArray<static_cast<std::size_t>(numComponents)>(
                fullyCompositionalThreePhase ? aqueousFugacityBegin : 0);

        // Preferred full-model naming.  `aqueousFugacity` remains an alias.
        static constexpr int waterFugacityBegin = aqueousFugacityBegin;
        static constexpr auto waterFugacity = aqueousFugacity;

        static constexpr int volumeClosure =
            fullyCompositionalThreePhase
                ? waterFugacityBegin + numComponents
                : fugacityBegin + numComponents;

        static constexpr int fullWellControl =
            fullyCompositionalThreePhase && hasWellUnknown
                ? volumeClosure + 1
                : disabledIndex;
        static constexpr int wellControl =
            fullyCompositionalThreePhase ? fullWellControl : legacyWellControl;

        static constexpr int legacyBaseCount = volumeClosure + 1;
        static constexpr int legacyAqueousCO2Equilibrium =
            !fullyCompositionalThreePhase && hasAqueousCO2Dissolution
                ? legacyBaseCount
                : disabledIndex;
        static constexpr int aqueousCO2Equilibrium = legacyAqueousCO2Equilibrium;

        static constexpr int count = fullyCompositionalThreePhase
            ? volumeClosure + 1 + (hasWellUnknown ? 1 : 0)
            : legacyBaseCount + (hasAqueousCO2Dissolution ? 1 : 0);
    };

    struct PhaseState final
    {
        static constexpr int flag = 0;

        // Legacy K_v/l block; in full mode this is also K_g/o.
        static constexpr int equilibriumRatioBegin = flag + 1;
        static constexpr auto equilibriumRatio =
            detail::makeIndexArray<static_cast<std::size_t>(numComponents)>(equilibriumRatioBegin);
        static constexpr auto vaporOilEquilibriumRatio = equilibriumRatio;

        static constexpr int aqueousOilEquilibriumRatioBegin =
            fullyCompositionalThreePhase
                ? equilibriumRatioBegin + numComponents
                : disabledIndex;
        static constexpr auto aqueousOilEquilibriumRatio =
            detail::makeIndexArray<static_cast<std::size_t>(numComponents)>(
                fullyCompositionalThreePhase ? aqueousOilEquilibriumRatioBegin : 0);
        static constexpr int waterOilEquilibriumRatioBegin = aqueousOilEquilibriumRatioBegin;
        static constexpr auto waterOilEquilibriumRatio = aqueousOilEquilibriumRatio;

        static constexpr int overallCompositionBegin =
            fullyCompositionalThreePhase
                ? aqueousOilEquilibriumRatioBegin + numComponents
                : equilibriumRatioBegin + numComponents;
        static constexpr auto overallComposition =
            detail::makeIndexArray<static_cast<std::size_t>(numComponents)>(overallCompositionBegin);

        // Legacy scalar L location remains unchanged in legacy mode.
        static constexpr int liquidFraction = overallCompositionBegin + numComponents;

        static constexpr int phaseMoleFractionBegin =
            fullyCompositionalThreePhase ? liquidFraction : disabledIndex;
        static constexpr auto phaseMoleFraction =
            detail::makeIndexArray<3>(fullyCompositionalThreePhase ? phaseMoleFractionBegin : 0);

        static constexpr int liquidCompressibilityFactor =
            fullyCompositionalThreePhase
                ? phaseMoleFractionBegin + 3
                : liquidFraction + 1;
        static constexpr int vaporCompressibilityFactor = liquidCompressibilityFactor + 1;
        static constexpr int waterCompressibilityFactor =
            fullyCompositionalThreePhase ? vaporCompressibilityFactor + 1 : disabledIndex;

        // 仅全组分三相使用：记录由微量饱和度 active-set 主动删除的相。放在
        // phase-state 块末尾，保持此前所有公开列/索引位置不变。
        static constexpr int phaseSuppressionFlag =
            fullyCompositionalThreePhase ? waterCompressibilityFactor + 1 : disabledIndex;

        static constexpr int count = fullyCompositionalThreePhase
            ? phaseSuppressionFlag + 1
            : vaporCompressibilityFactor + 1;
    };

    static constexpr int numPrimaryVariables = Primary::count;
    static constexpr int numEquations = Equation::count;
    static constexpr int numPhaseStateVariables = PhaseState::count;

    static_assert(
        numPrimaryVariables == numEquations,
        "Natural formulation must remain a square nonlinear system.");

    using ValueType = std::conditional_t<
        useAutomaticDifferentiation,
        DenseAd::Evaluation<double, numPrimaryVariables>,
        double>;
};

template <class Config>
using ADIndices = Indices<Config, true>;

template <class Config>
using ScalarIndices = Indices<Config, false>;

} // namespace MPMC
