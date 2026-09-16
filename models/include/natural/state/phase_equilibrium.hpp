/**
 * @file phase_equilibrium.hpp
 * @brief 两相稳定性与相平衡状态更新算法。
 */
#pragma once

#include <common/units.hpp>
#include <natural/fluid_system.hpp>
#include <natural/numerics.hpp>
#include <natural/phase_state.hpp>
#include <natural/state/phase_state_data.hpp>
#include <natural/state/phase_update_result.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace MPMC
{

namespace phase_detail
{
template <std::size_t N>
void normalize(std::array<double, N> &composition, double floor)
{
    double sum = 0.0;
    for (double &value : composition)
    {
        value = std::max(value, floor);
        sum += value;
    }
    if (!(sum > 0.0) || !std::isfinite(sum))
        throw std::runtime_error("Phase composition normalization failed.");
    for (double &value : composition)
        value /= sum;
}

template <class Indices>
std::array<double, Indices::numComponents>
compositionFromPrimary(
    const std::array<double, Indices::numPrimaryVariables> &primary,
    bool liquid)
{
    std::array<double, Indices::numComponents> result{};
    result.back() = 1.0;
    const auto &indices = liquid
        ? Indices::Primary::liquidComposition
        : Indices::Primary::vaporComposition;
    for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
    {
        result[static_cast<std::size_t>(i)] = primary[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])];
        result.back() -= result[static_cast<std::size_t>(i)];
    }
    return result;
}

template <class Indices>
void writeCompositionToPrimary(
    std::array<double, Indices::numPrimaryVariables> &primary,
    const std::array<double, Indices::numComponents> &composition,
    bool liquid)
{
    const auto &indices = liquid
        ? Indices::Primary::liquidComposition
        : Indices::Primary::vaporComposition;
    for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
        primary[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])] = composition[static_cast<std::size_t>(i)];
}
} // namespace phase_detail

/**
 * @brief 烃类相稳定测试、相出现/消失与 secondary phase-state 更新。
 *
 * Michelsen-style trial 以总体组成 `z_i` 和平衡比 `K_i=y_i/x_i` 为基础：
 * vapor-like trial 使用 `w_i=K_i z_i`，liquid-like trial 使用 `w_i=z_i/K_i`，
 * 归一化后通过两相逸度比做 successive substitution。
 *
 * 两相状态满足总体组成关系 `z_i=L x_i+(1-L)y_i`。相消失后，本类保持
 * 方程数不变并更新 K、z、液相摩尔比例 L 和两相 Z 因子，供下一 Newton
 * 线性化和下一时间步的稳定性测试使用。
 */
template <class Indices>
class PhaseEquilibriumManager final
{
public:
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;

    explicit PhaseEquilibriumManager(const FluidSystem<Indices> &fluid)
        : fluid_(fluid)
    {
    }

    /** @brief 一次气相型或液相型 Michelsen 稳定性试探结果。 */
    struct StabilityTrial
    {
        Composition trialComposition{};
        Composition equilibriumRatio{};
        double sum{0.0};
        bool trivial{false};
        bool converged{false};
        bool valid{false};
    };

    /**
     * @brief 执行一次完整的 Michelsen 型稳定性试探。
     *
     * `sum` 为未归一化试探组成之和；若非平凡试探收敛且 `sum>1`，说明该分相方向可降低 Gibbs 能。
     */
    [[nodiscard]] StabilityTrial stabilityTrial(
        const Composition &zInput,
        const Composition &initialK,
        double pressure,
        bool vaporLikeTrial) const
    {
        Composition z = zInput;
        phase_detail::normalize(z, NaturalNumerics::minimumComposition);
        Composition K = initialK;
        for (double &value : K)
            value = std::max(value, NaturalNumerics::minimumComposition);

        const auto feed = fluid_.eos.phaseResult(
            pressure, fluid_.temperature, z, vaporLikeTrial /* liquid root for vapor-like trial */);

        StabilityTrial result;
        for (int iteration = 0; iteration < NaturalNumerics::phaseStabilityMaximumIterations; ++iteration)
        {
            Composition trial{};
            double sum = 0.0;
            for (int component = 0; component < N; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                trial[c] = vaporLikeTrial ? K[c] * z[c] : z[c] / K[c];
                sum += trial[c];
            }
            if (!(sum > NaturalNumerics::minimumNormalizationDenominator) || !std::isfinite(sum))
                break;
            for (double &value : trial)
                value /= sum;

            const auto trialThermo = fluid_.eos.phaseResult(
                pressure, fluid_.temperature, trial, !vaporLikeTrial);

            double correctionNorm = 0.0;
            double kNorm = 0.0;
            for (int component = 0; component < N; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                double R;
                if (vaporLikeTrial)
                    R = feed.fugacity[c] / (trialThermo.fugacity[c] * sum);
                else
                    R = trialThermo.fugacity[c] / feed.fugacity[c] * sum;
                if (!(R > 0.0) || !std::isfinite(R))
                    return result;
                K[c] *= R;
                correctionNorm += (R - 1.0) * (R - 1.0);
                kNorm += std::log(K[c]) * std::log(K[c]);
            }

            result.trialComposition = trial;
            result.equilibriumRatio = K;
            result.sum = sum;
            result.trivial = kNorm < 1.0e-5;
            result.converged = correctionNorm < NaturalNumerics::phaseStabilityTolerance;
            if (result.trivial || result.converged)
            {
                result.valid = true;
                return result;
            }
        }
        return result;
    }

    /** @brief 汇总液相型/气相型试探后的稳定性结论及候选组成。 */
    struct StabilityResult
    {
        bool stable{true};
        Composition liquid{};
        Composition vapor{};
    };

    /** @brief 判断当前烃类总体组成作为单相是否热力学稳定。 */
    [[nodiscard]] StabilityResult testStability(
        Composition z,
        const Composition &K,
        double pressure) const
    {
        phase_detail::normalize(z, NaturalNumerics::minimumComposition);

        if (fluid_.eos.usesEquilibriumConstants())
        {
            const auto updatedK = fluid_.eos.evaluateEquilibriumConstants(
                pressure, fluid_.temperature, z);
            const double liquidFraction = fluid_.eos.solveRachfordRice(
                0.0, updatedK, z);
            constexpr double tolerance = 1.0e-10;
            const bool stable =
                std::abs(liquidFraction - 1.0) <= tolerance || liquidFraction <= tolerance;
            return {
                stable,
                fluid_.eos.liquidComposition(liquidFraction, updatedK, z),
                fluid_.eos.vaporComposition(liquidFraction, updatedK, z)};
        }

        const auto vaporTrial = stabilityTrial(z, K, pressure, true);
        const auto liquidTrial = stabilityTrial(z, K, pressure, false);

        /*
         * A failed local stability iteration is not a thermodynamic statement.
         * In particular, StabilityTrial::sum defaults to zero, so using only
         * `sum <= 1` would accidentally classify a failed trial as stable.
         * Keep the current single phase conservatively when either trial is
         * invalid; the global Newton state remains untouched and a later state
         * update can retry with a different p/z/K.
         */
        if (!vaporTrial.valid || !liquidTrial.valid)
            return {true, z, z};

        // Pre-activate a non-trivial incipient phase inside the same 1e-4
        // hysteresis band used by the fully-compositional active set.  Waiting
        // until sum is already above one makes the legacy variable switch
        // coincide with a sharp transport front and produces a non-smooth
        // Newton branch change.  Trivial trials remain single-phase.
        constexpr double appearanceThreshold =
            1.0 - NaturalNumerics::phaseHysteresisReappearanceMargin;
        const bool vaporStable =
            vaporTrial.trivial || vaporTrial.sum <= appearanceThreshold;
        const bool liquidStable =
            liquidTrial.trivial || liquidTrial.sum <= appearanceThreshold;
        StabilityResult result;
        result.stable = vaporStable && liquidStable;
        if (result.stable)
        {
            result.liquid = z;
            result.vapor = z;
        }
        else
        {
            result.liquid = liquidTrial.trialComposition;
            result.vapor = vaporTrial.trialComposition;
            if (!std::isfinite(liquidTrial.sum)) result.liquid = z;
            if (!std::isfinite(vaporTrial.sum)) result.vapor = z;
        }
        return result;
    }

    /**
     * @brief 更新 `K_i`、`z_i`、烃液相摩尔分率 `L` 和两相 Z 因子。
     *
     * From phase saturation and molar density,
     * `L = rho_m,l S_l / (rho_m,l S_l + rho_m,g S_g)`. Then
     * `z_i=L x_i+(1-L)y_i` and, where `x_i>0`, `K_i=y_i/x_i`.
     */
    void updateSecondary(
        std::array<double, Indices::numPrimaryVariables> &primary,
        PhaseStateData<Indices> &phaseState) const
    {
        auto x = phase_detail::compositionFromPrimary<Indices>(primary, true);
        auto y = phase_detail::compositionFromPrimary<Indices>(primary, false);
        phase_detail::normalize(x, NaturalNumerics::minimumComposition);
        phase_detail::normalize(y, NaturalNumerics::minimumComposition);

        double sL = primary[Indices::Primary::liquidSaturation];
        double sV = primary[Indices::Primary::vaporSaturation];
        const double hydrocarbonSaturation = sL + sV;
        if (hydrocarbonSaturation > NaturalNumerics::minimumNormalizationDenominator)
        {
            sL /= hydrocarbonSaturation;
            sV /= hydrocarbonSaturation;
        }
        else
        {
            sL = std::clamp(phaseState.liquidMoleFraction, 0.0, 1.0);
            sV = 1.0 - sL;
        }

        const double p = primary[Indices::Primary::pressure];
        double rhoMolL = 0.0;
        double rhoMolV = 0.0;
        if (!fluid_.eos.usesEquilibriumConstants())
        {
            const auto l = fluid_.eos.phaseResult(p, fluid_.temperature, x, true);
            const auto v = fluid_.eos.phaseResult(p, fluid_.temperature, y, false);
            phaseState.liquidCompressibility = l.compressibility;
            phaseState.vaporCompressibility = v.compressibility;
            // Keep secondary-state molar density on the exact EOS path.  The
            // helper handles both translated and non-translated cubic models and
            // therefore also guarantees a single gas-constant convention.
            rhoMolL = fluid_.eos.molarDensity(
                p, fluid_.temperature, x, l.compressibility,
                CompositionalPhase::Oil);
            rhoMolV = fluid_.eos.molarDensity(
                p, fluid_.temperature, y, v.compressibility,
                CompositionalPhase::Gas);
        }
        else
        {
            rhoMolL = fluid_.blackOilProperties.computeMolarDensity(p, x, true);
            rhoMolV = fluid_.blackOilProperties.computeMolarDensity(p, y, false);
            phaseState.liquidCompressibility =
                p / (units::gasConstant * fluid_.temperature * rhoMolL);
            phaseState.vaporCompressibility =
                p / (units::gasConstant * fluid_.temperature * rhoMolV);
        }
        const double denominator = rhoMolL * sL + rhoMolV * sV;
        if (denominator > NaturalNumerics::minimumNormalizationDenominator)
            phaseState.liquidMoleFraction = rhoMolL * sL / denominator;

        /*
         * 数值：`z_i` 在单/两相状态都有效，但 `K_i=y_i/x_i` 只在两相真实共存时有物理意义。
         * 单相时缺失相组成只是有限占位值；若此时重算 K，会把所有 K 强制成 1，导致下一次
         * Michelsen 稳定性测试从平凡解出发而无法发现新相。因此单相阶段保留最后一次物理有效 K。
         */
        const bool updateEquilibriumRatio =
            phaseState.phase == HydrocarbonPhaseState::TwoPhase;

        for (int component = 0; component < N; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            phaseState.overallComposition[c] =
                phaseState.liquidMoleFraction * x[c] +
                (1.0 - phaseState.liquidMoleFraction) * y[c];
            if (updateEquilibriumRatio &&
                x[c] > NaturalNumerics::minimumNormalizationDenominator)
            {
                phaseState.equilibriumRatio[c] = y[c] / x[c];
            }
        }
    }

    /**
     * @brief 用稳定性测试和饱和度判据处理烃相出现/消失。
     *
     * 状态：单相只有在稳定性测试判为不稳定后才生成缺失相；两相中饱和度降到消失阈值的相被移除，
     * 其物质通过总体组成继续保留在剩余相平衡中。
     */
    PhaseUpdateResult updatePhaseState(
        std::array<double, Indices::numPrimaryVariables> &primary,
        PhaseStateData<Indices> &phaseState) const
    {
        auto x = phase_detail::compositionFromPrimary<Indices>(primary, true);
        auto y = phase_detail::compositionFromPrimary<Indices>(primary, false);
        phase_detail::normalize(x, NaturalNumerics::minimumComposition);
        phase_detail::normalize(y, NaturalNumerics::minimumComposition);

        Composition z{};
        double zSum = 0.0;
        for (int component = 0; component < N; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            if (phaseState.phase == HydrocarbonPhaseState::LiquidOnly)
                z[c] = x[c];
            else if (phaseState.phase == HydrocarbonPhaseState::VaporOnly)
                z[c] = y[c];
            else
                z[c] = phaseState.liquidMoleFraction * x[c] +
                       (1.0 - phaseState.liquidMoleFraction) * y[c];
            zSum += z[c];
        }
        if (!(zSum > 0.0))
            throw std::runtime_error("Overall composition sum must be positive.");
        for (double &value : z) value /= zSum;
        phaseState.overallComposition = z;

        bool stable = phaseState.phase != HydrocarbonPhaseState::TwoPhase;
        if (phaseState.phase != HydrocarbonPhaseState::TwoPhase)
        {
            const auto stability = testStability(
                z, phaseState.equilibriumRatio, primary[Indices::Primary::pressure]);
            stable = stability.stable;
            x = stability.liquid;
            y = stability.vapor;
        }

        double sL = primary[Indices::Primary::liquidSaturation];
        double sV = primary[Indices::Primary::vaporSaturation];
        const bool toEpsilonLiquid =
            phaseState.phase == HydrocarbonPhaseState::VaporOnly && !stable;
        const bool toEpsilonVapor =
            phaseState.phase == HydrocarbonPhaseState::LiquidOnly && !stable;
        const bool toOnlyLiquid =
            phaseState.phase == HydrocarbonPhaseState::TwoPhase && sV <= 0.0;
        const bool toOnlyVapor =
            phaseState.phase == HydrocarbonPhaseState::TwoPhase && sL <= 0.0;

        // A stability test returns equilibrium-like compositions for both
        // trial phases.  At phase appearance, however, only the incipient
        // phase may take its trial composition: replacing the existing
        // single-phase composition would create an O(1) inventory jump while
        // the new phase carries only the seed saturation.  Retain z in the
        // existing phase and let Newton move both compositions continuously.
        if (toEpsilonLiquid)
            y = z;
        if (toEpsilonVapor)
            x = z;

        const bool pureLiquid =
            (stable && phaseState.phase == HydrocarbonPhaseState::LiquidOnly) || toOnlyLiquid;
        const bool pureVapor =
            (stable && phaseState.phase == HydrocarbonPhaseState::VaporOnly && !pureLiquid) || toOnlyVapor;

        if (toEpsilonLiquid || toEpsilonVapor)
        {
            for (int component = 0; component < N; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                phaseState.equilibriumRatio[c] = y[c] / x[c];
            }
        }

        double nonHydrocarbonSaturation = 0.0;
        if constexpr (Indices::hasWater)
            nonHydrocarbonSaturation = primary[Indices::Primary::waterSaturation];
        const double maximumHydrocarbonSaturation =
            std::max(1.0 - nonHydrocarbonSaturation, NaturalNumerics::minimumComposition);
        const double seed = NaturalNumerics::phaseAppearanceSaturation;

        if (toEpsilonLiquid)
        {
            sL = seed;
            sV = maximumHydrocarbonSaturation - seed;
        }
        if (toEpsilonVapor)
        {
            sV = seed;
            sL = maximumHydrocarbonSaturation - seed;
        }
        if (pureVapor)
        {
            sL = 0.0;
            sV = maximumHydrocarbonSaturation;
            x = y;
            phaseState.phase = HydrocarbonPhaseState::VaporOnly;
        }
        else if (pureLiquid)
        {
            sL = maximumHydrocarbonSaturation;
            sV = 0.0;
            y = x;
            phaseState.phase = HydrocarbonPhaseState::LiquidOnly;
        }
        else
        {
            phaseState.phase = HydrocarbonPhaseState::TwoPhase;
        }

        primary[Indices::Primary::liquidSaturation] = sL;
        primary[Indices::Primary::vaporSaturation] = sV;
        phase_detail::normalize(x, NaturalNumerics::minimumComposition);
        phase_detail::normalize(y, NaturalNumerics::minimumComposition);
        phase_detail::writeCompositionToPrimary<Indices>(primary, x, true);
        phase_detail::writeCompositionToPrimary<Indices>(primary, y, false);
        updateSecondary(primary, phaseState);
        return {
            phaseState.phase == HydrocarbonPhaseState::TwoPhase
                ? PhaseUpdateStatus::Updated
                : PhaseUpdateStatus::StableReducedSet};
    }

    /**
     * @brief updateState 回调中最前端的状态清理/归一化。
     */
    void sanitizePrimaryBeforeFlash(
        std::array<double, Indices::numPrimaryVariables> &primary,
        bool useVariableBounds) const
    {
        if (!useVariableBounds)
        {
            primary[Indices::Primary::liquidSaturation] =
                std::clamp(primary[Indices::Primary::liquidSaturation], 0.0, 1.0);
            primary[Indices::Primary::vaporSaturation] =
                std::clamp(primary[Indices::Primary::vaporSaturation], 0.0, 1.0);
            if constexpr (Indices::hasWater)
                primary[Indices::Primary::waterSaturation] =
                    std::clamp(primary[Indices::Primary::waterSaturation], 0.0, 1.0);
            if constexpr (Indices::hasAqueousCO2Dissolution)
                primary[Indices::Primary::aqueousCO2MoleFraction] =
                    std::clamp(primary[Indices::Primary::aqueousCO2MoleFraction], 0.0, 1.0 - 1.0e-12);
        }

        // 这里只裁剪物理边界，不额外归一化饱和度；归一化由相平衡/闭合方程
        // 自身完成，避免预处理改变 Newton 线性化。
    }

private:
    const FluidSystem<Indices> &fluid_;
};

} // namespace MPMC
