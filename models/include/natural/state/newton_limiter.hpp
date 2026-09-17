/**
 * @file newton_limiter.hpp
 * @brief Newton 更新的压力、饱和度和组成步长限制。
 */
#pragma once

#include <natural/numerics.hpp>
#include <natural/phase_state.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace MPMC
{

namespace detail
{
inline double limiterScale(double magnitude, double maximum) noexcept
{
    if (!(magnitude > 0.0) || !std::isfinite(magnitude))
        return 1.0;
    return std::min(maximum / magnitude, 1.0);
}
} // namespace detail

/**
 * @brief 限制单元 Newton 增量，并显式处理零尺度和非有限值。
 *
 * pressure 单独按相对变化限制；饱和度、两相独立组成和 aqueous CO2 共用
 * thermodynamic scale，使各热力学未知量同步缩放。
 */
template <class Indices>
void limitNaturalNewtonIncrement(
    std::array<double, Indices::numPrimaryVariables> &delta,
    const std::array<double, Indices::numPrimaryVariables> &state,
    HydrocarbonPhaseState phaseState)
{
    if (phaseState == HydrocarbonPhaseState::VaporOnly)
    {
        for (int component = 0; component < Indices::numIndependentCompositionsPerPhase; ++component)
        {
            const int liq = Indices::Primary::liquidComposition[static_cast<std::size_t>(component)];
            const int vap = Indices::Primary::vaporComposition[static_cast<std::size_t>(component)];
            delta[static_cast<std::size_t>(vap)] = delta[static_cast<std::size_t>(liq)];
            delta[static_cast<std::size_t>(liq)] = 0.0;
        }
    }

    double maximumSaturationDelta =
        std::max(std::abs(delta[Indices::Primary::liquidSaturation]),
                 std::abs(delta[Indices::Primary::vaporSaturation]));
    if constexpr (Indices::hasWater)
        maximumSaturationDelta = std::max(
            maximumSaturationDelta,
            std::abs(delta[Indices::Primary::waterSaturation]));
    const double saturationScale = detail::limiterScale(
        maximumSaturationDelta, NaturalNumerics::maximumSaturationNewtonChange);

    double maxLiquidCompositionDelta = 0.0;
    double maxVaporCompositionDelta = 0.0;
    double dependentLiquidDelta = 0.0;
    double dependentVaporDelta = 0.0;
    for (int component = 0; component < Indices::numIndependentCompositionsPerPhase; ++component)
    {
        const double dx = delta[Indices::Primary::liquidComposition[static_cast<std::size_t>(component)]];
        const double dy = delta[Indices::Primary::vaporComposition[static_cast<std::size_t>(component)]];
        dependentLiquidDelta -= dx;
        dependentVaporDelta -= dy;
        maxLiquidCompositionDelta = std::max(maxLiquidCompositionDelta, std::abs(dx));
        maxVaporCompositionDelta = std::max(maxVaporCompositionDelta, std::abs(dy));
    }
    maxLiquidCompositionDelta = std::max(maxLiquidCompositionDelta, std::abs(dependentLiquidDelta));
    maxVaporCompositionDelta = std::max(maxVaporCompositionDelta, std::abs(dependentVaporDelta));
    const double compositionScale = std::min(
        detail::limiterScale(maxLiquidCompositionDelta, NaturalNumerics::maximumCompositionNewtonChange),
        detail::limiterScale(maxVaporCompositionDelta, NaturalNumerics::maximumCompositionNewtonChange));

    double aqueousScale = 1.0;
    if constexpr (Indices::hasAqueousCO2Dissolution)
    {
        aqueousScale = detail::limiterScale(
            std::abs(delta[Indices::Primary::aqueousCO2MoleFraction]),
            NaturalNumerics::maximumAqueousCO2NewtonChange);
    }

    const double thermodynamicScale =
        std::min({saturationScale, compositionScale, aqueousScale});
    delta[Indices::Primary::liquidSaturation] *= thermodynamicScale;
    delta[Indices::Primary::vaporSaturation] *= thermodynamicScale;
    if constexpr (Indices::hasWater)
        delta[Indices::Primary::waterSaturation] *= thermodynamicScale;
    for (int component = 0; component < Indices::numIndependentCompositionsPerPhase; ++component)
    {
        delta[Indices::Primary::liquidComposition[static_cast<std::size_t>(component)]] *= thermodynamicScale;
        delta[Indices::Primary::vaporComposition[static_cast<std::size_t>(component)]] *= thermodynamicScale;
    }
    if constexpr (Indices::hasAqueousCO2Dissolution)
        delta[Indices::Primary::aqueousCO2MoleFraction] *= thermodynamicScale;

    const double pressure = std::abs(state[Indices::Primary::pressure]);
    if (pressure > NaturalNumerics::minimumNormalizationDenominator)
    {
        const double relativeChange =
            std::abs(delta[Indices::Primary::pressure]) / pressure;
        delta[Indices::Primary::pressure] *= detail::limiterScale(
            relativeChange, NaturalNumerics::maximumRelativePressureNewtonChange);
    }
}

/**
 * @brief 全组分三相模式的 phase-aware Newton 增量限制器。
 *
 * 三个饱和度仍共享一个缩放因子，以保持 saturation-closure 的 Newton 方向；
 * 但 O/G/W 三个组成块分别限制。此前六个热力学块共用全局最小 scale，导致
 * disappearing phase 的病态组成增量可以把另外两相和全部饱和度一起压到近零，
 * 在相边界形成“线性求解收敛、非线性残差完全不动”的固定点。
 *
 * 组成块内部仍保持统一 scale，因而 N-1 独立摩尔分数与依赖末组分的方向一致；
 * 不同相之间无需共享这个人为的 damping。饱和度也故意不截断为非负，原始
 * Newton 信号交给 FullyCompositionalThreePhaseEquilibrium::updatePhaseState()
 * 统一处理相出现/消失。
 */
template <class Indices>
void limitNaturalNewtonIncrement(
    std::array<double, Indices::numPrimaryVariables> &delta,
    const std::array<double, Indices::numPrimaryVariables> &state,
    PhasePresence /*phasePresence*/)
{
    static_assert(Indices::fullyCompositionalThreePhase,
                  "PhasePresence limiter is only for the full three-phase formulation.");

    // Saturations keep one common scale so the Newton direction of the volume
    // closure is preserved. Oil composition is stored as q_i=S_o x_i; apply the
    // same first scale to q so the phase-amount direction follows S_o, then use
    // only an additional q-only scale if the complete amount vector would move
    // by more than one saturation-sized Newton step. A large minority-phase q
    // correction therefore cannot freeze the three saturation updates.
    const double maximumSaturationDelta = std::max({
        std::abs(delta[Indices::Primary::liquidSaturation]),
        std::abs(delta[Indices::Primary::vaporSaturation]),
        std::abs(delta[Indices::Primary::waterSaturation])});
    const double saturationScale = detail::limiterScale(
        maximumSaturationDelta, NaturalNumerics::maximumSaturationNewtonChange);

    delta[Indices::Primary::liquidSaturation] *= saturationScale;
    delta[Indices::Primary::vaporSaturation] *= saturationScale;
    delta[Indices::Primary::waterSaturation] *= saturationScale;

    double maximumOilIndependentAmountDelta = 0.0;
    double oilIndependentAmountSumDelta = 0.0;
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const auto primary = static_cast<std::size_t>(
            Indices::Primary::liquidComposition[static_cast<std::size_t>(component)]);
        delta[primary] *= saturationScale;
        maximumOilIndependentAmountDelta = std::max(
            maximumOilIndependentAmountDelta, std::abs(delta[primary]));
        oilIndependentAmountSumDelta += delta[primary];
    }

    // q_N=S_o-sum(q_i).  Reserve enough of the 0.1 amount-step budget for the
    // already-limited dS_o so both independent q_i and dependent q_N remain
    // bounded without feeding the q pathology back into the saturation scale.
    const double oilSaturationDelta =
        delta[Indices::Primary::liquidSaturation];
    const double remainingOilAmountBudget = std::max(
        0.0, NaturalNumerics::maximumSaturationNewtonChange -
                 std::abs(oilSaturationDelta));
    const double oilAmountScale = detail::limiterScale(
        std::max(maximumOilIndependentAmountDelta,
                 std::abs(oilIndependentAmountSumDelta)),
        remainingOilAmountBudget);
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        delta[static_cast<std::size_t>(
            Indices::Primary::liquidComposition[static_cast<std::size_t>(component)])] *=
            oilAmountScale;
    }

    const auto limitCompositionBlock = [&](const auto &indices)
    {
        double maximum = 0.0;
        double dependent = 0.0;
        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            const double value = delta[static_cast<std::size_t>(
                indices[static_cast<std::size_t>(component)])];
            maximum = std::max(maximum, std::abs(value));
            dependent -= value;
        }
        maximum = std::max(maximum, std::abs(dependent));
        const double scale = detail::limiterScale(
            maximum, NaturalNumerics::maximumCompositionNewtonChange);

        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            delta[static_cast<std::size_t>(
                indices[static_cast<std::size_t>(component)])] *= scale;
        }
    };

    // Gas and water stay ordinary mole-fraction coordinates and keep their own
    // phase-local scales.
    limitCompositionBlock(Indices::Primary::vaporComposition);
    limitCompositionBlock(Indices::Primary::waterComposition);

    const double pressure = std::abs(state[Indices::Primary::pressure]);
    if (pressure > NaturalNumerics::minimumNormalizationDenominator)
    {
        const double relativeChange =
            std::abs(delta[Indices::Primary::pressure]) / pressure;
        delta[Indices::Primary::pressure] *= detail::limiterScale(
            relativeChange, NaturalNumerics::maximumRelativePressureNewtonChange);
    }
}

/**
 * @brief 将全组分水相 Newton 候选态限制在 H2O/CO2 闭包的组分适用域内。
 *
 * 只缩放 water-composition 块并保持原 Newton 方向。末组分由
 * `x_N=1-sum(x_i)` 恢复，因此显式计入其当前值和增量。
 */
template <class Indices>
void limitAqueousCompositionNewtonIncrement(
    std::array<double, Indices::numPrimaryVariables> &delta,
    const std::array<double, Indices::numPrimaryVariables> &state,
    int waterComponent,
    int co2Component,
    double maximumUnsupportedMoleFraction)
{
    static_assert(Indices::fullyCompositionalThreePhase,
                  "Aqueous composition domain limiting requires the full three-phase formulation.");

    double unsupported = 0.0;
    double unsupportedDelta = 0.0;
    double dependent = 1.0;
    double dependentDelta = 0.0;

    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const auto primary = static_cast<std::size_t>(
            Indices::Primary::waterComposition[static_cast<std::size_t>(component)]);
        const double value = state[primary];
        const double change = delta[primary];
        dependent -= value;
        dependentDelta -= change;
        if (component != waterComponent && component != co2Component)
        {
            unsupported += value;
            unsupportedDelta += change;
        }
    }

    constexpr int dependentComponent = Indices::numComponents - 1;
    if (dependentComponent != waterComponent && dependentComponent != co2Component)
    {
        unsupported += dependent;
        unsupportedDelta += dependentDelta;
    }

    if (!(unsupportedDelta > 0.0))
        return;

    // Leave a tiny round-off margin so EOS/viscosity scalar checks never see a
    // value infinitesimally above the identical configured ceiling.
    const double margin = std::max(
        NaturalNumerics::phaseEquilibriumTraceComposition,
        maximumUnsupportedMoleFraction * 1.0e-10);
    const double target = std::max(0.0, maximumUnsupportedMoleFraction - margin);
    const double scale = std::clamp(
        (target - unsupported) / unsupportedDelta, 0.0, 1.0);

    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        delta[static_cast<std::size_t>(
            Indices::Primary::waterComposition[static_cast<std::size_t>(component)])] *= scale;
    }
}

} // namespace MPMC
