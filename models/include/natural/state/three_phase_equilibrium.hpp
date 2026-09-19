/**
 * @file three_phase_equilibrium.hpp
 * @brief O/G/W 三相状态的相稳定性判断与 flash 调度。
 */
#pragma once

#include <natural/fluid_system.hpp>
#include <natural/numerics.hpp>
#include <natural/state/phase_state_data.hpp>
#include <natural/state/phase_update_result.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

namespace three_phase_detail
{
template <std::size_t N>
void normalize(std::array<double, N> &composition, double floor)
{
    double sum = 0.0;
    for (double &v : composition)
    {
        v = std::max(v, floor);
        sum += v;
    }
    if (!(sum > 0.0) || !std::isfinite(sum))
        throw std::runtime_error("Fully compositional phase composition normalization failed.");
    for (double &v : composition) v /= sum;
}

template <class Indices>
std::array<double, Indices::numComponents> compositionFromPrimary(
    const std::array<double, Indices::numPrimaryVariables> &primary,
    CompositionalPhase phase)
{
    std::array<double, Indices::numComponents> result{};
    const auto *indices = &Indices::Primary::liquidComposition;
    if (phase == CompositionalPhase::Gas)
        indices = &Indices::Primary::vaporComposition;
    else if (phase == CompositionalPhase::Water)
        indices = &Indices::Primary::waterComposition;

    result.back() = 1.0;
    for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
    {
        const std::size_t c = static_cast<std::size_t>(i);
        result[c] = primary[static_cast<std::size_t>((*indices)[c])];
        result.back() -= result[c];
    }
    return result;
}

template <class Indices>
void writeCompositionToPrimary(
    std::array<double, Indices::numPrimaryVariables> &primary,
    const std::array<double, Indices::numComponents> &composition,
    CompositionalPhase phase)
{
    const auto *indices = &Indices::Primary::liquidComposition;
    if (phase == CompositionalPhase::Gas)
        indices = &Indices::Primary::vaporComposition;
    else if (phase == CompositionalPhase::Water)
        indices = &Indices::Primary::waterComposition;

    for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
    {
        const std::size_t c = static_cast<std::size_t>(i);
        primary[static_cast<std::size_t>((*indices)[c])] = composition[c];
    }
}

inline int saturationIndex(CompositionalPhase phase, int liquid, int vapor, int water)
{
    switch (phase)
    {
    case CompositionalPhase::Oil: return liquid;
    case CompositionalPhase::Gas: return vapor;
    case CompositionalPhase::Water: return water;
    }
    return liquid;
}
} // namespace three_phase_detail

/**
 * @brief 全组分 O/G/W 相态管理器，连接 Newton 主变量与热力学 flash。
 *
 * 所有 N 个组分（包括 H2O）均允许进入三相。主变量保存三相组成和饱和度。
 * 状态：Newton 状态若把某相饱和度推进到 active-set 边界容差，该相先被约化移除，
 * 剩余相集再执行与 P–T–z flash 相同的多相 TPD 稳定性测试；缺失相仅在越过
 * 稳定性滞回带后重新生成，避免相边界附近的 active-set 抖动。
 */
template <class Indices, class FlashBackend = CubicThreePhaseFlash<Indices>>
class FullyCompositionalThreePhaseEquilibrium final
{
public:
    static_assert(Indices::fullyCompositionalThreePhase,
                  "FullyCompositionalThreePhaseEquilibrium requires the full three-phase model config.");
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;
    using PrimaryArray = std::array<double, Indices::numPrimaryVariables>;
    using Flash = FlashBackend;
    using FlashResult = typename Flash::Result;
    using StabilityResult = typename Flash::StabilityResult;

    explicit FullyCompositionalThreePhaseEquilibrium(const FluidSystem<Indices> &fluid)
        : fluid_(fluid),
          flash_(fluid.eos, flashOptions_(fluid))
    {
        if (fluid.eos.usesEquilibriumConstants())
            throw std::logic_error("Fully compositional three-phase mode requires an EOS fugacity model, not K-value tables.");
    }

    [[nodiscard]] FlashResult flashPTZ(
        double pressure,
        double temperature,
        Composition overallComposition) const
    {
        return flash_.flash(pressure, temperature, std::move(overallComposition));
    }

    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */
    void assignFlashResult(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState,
        const FlashResult &flash) const
    {
        if (!flash.converged)
            throw std::runtime_error("Fully compositional PTz flash did not converge.");

        primary[Indices::Primary::liquidSaturation] = flash.saturation[0];
        primary[Indices::Primary::vaporSaturation] = flash.saturation[1];
        primary[Indices::Primary::waterSaturation] = flash.saturation[2];
        for (CompositionalPhase phase : phases_)
            three_phase_detail::writeCompositionToPrimary<Indices>(
                primary,
                flash.composition[static_cast<std::size_t>(phaseIndex(phase))],
                phase);

        phaseState.phasePresence = flash.presence;
        for (CompositionalPhase phase : phases_)
        {
            if (flash.presence.contains(phase))
                phaseState.phaseSuppression.remove(phase);
        }
        phaseState.overallComposition = overallFromFlash_(flash);
        phaseState.vaporOilEquilibriumRatio = flash.vaporOilK;
        phaseState.waterOilEquilibriumRatio = flash.waterOilK;
        phaseState.phaseMoleFraction = flash.phaseMoleFraction;
        phaseState.phaseCompressibility = flash.compressibility;
    }

    /**
     * @brief 由当前主变量更新 secondary `z`、`beta`、K 和 Z。
     *
     * 数学：`beta_a = rho_m,a S_a / sum(rho_m S)`。未活动相的 K 值保留最后一次物理有效值，
     * 避免单相占位组成把 K 强制变成 1，从而使后续稳定性测试锁死在平凡解。
     */
    void updateSecondary(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState) const
    {
        const auto composition = normalizedPhaseCompositions_(primary);
        std::array<double, 3> saturation{
            primary[Indices::Primary::liquidSaturation],
            primary[Indices::Primary::vaporSaturation],
            primary[Indices::Primary::waterSaturation]};
        std::array<double, 3> molarDensity{};
        std::array<double, 3> moleAmount{};

        const double p = primary[Indices::Primary::pressure];
        double totalMoles = 0.0;

        const auto accumulatePhaseMoles = [&](CompositionalPhase phase, double compressibility) {
            const std::size_t pi = static_cast<std::size_t>(phaseIndex(phase));
            phaseState.phaseCompressibility[pi] = compressibility;
            molarDensity[pi] = fluid_.eos.molarDensity(
                p, fluid_.temperature, composition[pi], compressibility, phase);
            moleAmount[pi] = molarDensity[pi] * std::max(saturation[pi], 0.0);
            totalMoles += moleAmount[pi];
        };

        if (fluid_.eos.usesCubicPlusAssociation())
        {
            for (CompositionalPhase phase : phases_)
            {
                if (!phaseState.phasePresence.contains(phase))
                    continue;
                const std::size_t pi = static_cast<std::size_t>(phaseIndex(phase));
                accumulatePhaseMoles(
                    phase,
                    fluid_.eos.compressibility(
                        p, fluid_.temperature, composition[pi], phase));
            }
        }
        else
        {
            // 性能：普通 PR 各相共享一套 (p,T) 参数；SW 的油/气共享非水参数，只额外需要水相参数。
            // 复用参数矩阵只消除重复 alpha/BIP 计算，不改变组成混合和立方根求解。
            const auto nonAqueousParameters = fluid_.eos.mixingParameters(
                p, fluid_.temperature, CompositionalPhase::Oil);
            const auto aqueousParameters = fluid_.eos.usesSoreideWhitson()
                ? fluid_.eos.mixingParameters(
                      p, fluid_.temperature, CompositionalPhase::Water)
                : nonAqueousParameters;

            for (CompositionalPhase phase : phases_)
            {
                if (!phaseState.phasePresence.contains(phase))
                    continue;
                const std::size_t pi = static_cast<std::size_t>(phaseIndex(phase));
                const auto &parameters = phase == CompositionalPhase::Water
                    ? aqueousParameters
                    : nonAqueousParameters;
                const auto mix = fluid_.eos.phaseMixing(composition[pi], parameters);
                const double z = phase == CompositionalPhase::Gas
                    ? fluid_.eos.vaporRoot(mix.A, mix.B)
                    : fluid_.eos.liquidRoot(mix.A, mix.B);
                accumulatePhaseMoles(phase, z);
            }
        }

        if (!(totalMoles > NaturalNumerics::minimumNormalizationDenominator))
            throw std::runtime_error("Fully compositional secondary update has zero active phase moles.");

        phaseState.overallComposition.fill(0.0);
        for (int pidx = 0; pidx < 3; ++pidx)
        {
            phaseState.phaseMoleFraction[static_cast<std::size_t>(pidx)] =
                moleAmount[static_cast<std::size_t>(pidx)] / totalMoles;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                phaseState.overallComposition[c] +=
                    phaseState.phaseMoleFraction[static_cast<std::size_t>(pidx)] *
                    composition[static_cast<std::size_t>(pidx)][c];
            }
        }
        three_phase_detail::normalize(
            phaseState.overallComposition,
            flash_.options().compositionFloor);

        if (phaseState.phasePresence.contains(CompositionalPhase::Oil) &&
            phaseState.phasePresence.contains(CompositionalPhase::Gas))
        {
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                phaseState.vaporOilEquilibriumRatio[c] =
                    composition[1][c] /
                    std::max(composition[0][c], flash_.options().compositionFloor);
            }
        }
        if (phaseState.phasePresence.contains(CompositionalPhase::Oil) &&
            phaseState.phasePresence.contains(CompositionalPhase::Water))
        {
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                phaseState.waterOilEquilibriumRatio[c] =
                    composition[2][c] /
                    std::max(composition[0][c], flash_.options().compositionFloor);
            }
        }
    }

    /**
     * @brief 先移除低于 active-set 容差的相，再用带滞回的稳定性测试判断是否重新生成。
     *
     * The overall composition used for the phase transition is captured before
     * any disappearing-phase saturation is reset.  This is essential: zeroing a phase
     * first and then rebuilding z would silently delete the material carried by
     * the disappearing phase.  We instead combine the current Newton
     * compositions with the last physical phase mole fractions, mirroring the
     * variable-switching logic of the legacy oil/gas formulation.
     */
    PhaseUpdateResult updatePhaseState(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState) const
    {
        canonicalizeSwNonAqueousRoles_(primary, phaseState);
        const PrimaryArray originalPrimary = primary;
        const PhaseStateData<Indices> originalPhaseState = phaseState;
        PhaseUpdateResult updateResult{};

        const Composition transitionOverall =
            transitionOverallComposition_(primary, phaseState);

        PhasePresence active = phaseState.phasePresence;
        bool removed = false;
        const std::array<double, 3> newtonSaturation{
            primary[Indices::Primary::liquidSaturation],
            primary[Indices::Primary::vaporSaturation],
            primary[Indices::Primary::waterSaturation]};

        // phaseSuppression also carries the short-lived appearance hold for a
        // phase that stability has legitimately reintroduced inside the probe
        // band.  Do not clear that memory merely because the phase is active:
        // otherwise a flash returning S << phaseBoundaryProbeSaturation is
        // deleted again at the very next post-check and Newton chatters across
        // two different active sets.  Once the phase grows outside the probe
        // band it becomes an ordinary active phase and the hold is released.
        for (CompositionalPhase phase : phases_)
        {
            const std::size_t p =
                static_cast<std::size_t>(phaseIndex(phase));
            if (active.contains(phase) &&
                phaseState.phaseSuppression.contains(phase) &&
                newtonSaturation[p] >
                    NaturalNumerics::phaseBoundaryProbeSaturation)
            {
                phaseState.phaseSuppression.remove(phase);
            }
        }

        for (CompositionalPhase phase : phases_)
        {
            const int sIndex = three_phase_detail::saturationIndex(
                phase,
                Indices::Primary::liquidSaturation,
                Indices::Primary::vaporSaturation,
                Indices::Primary::waterSaturation);
            const double saturation =
                newtonSaturation[static_cast<std::size_t>(phaseIndex(phase))];
            const bool appearanceHold =
                phaseState.phaseSuppression.contains(phase) &&
                saturation > 0.0 &&
                saturation <= NaturalNumerics::phaseBoundaryProbeSaturation;
            if (active.contains(phase) &&
                saturation <= NaturalNumerics::phaseBoundaryProbeSaturation &&
                !appearanceHold)
            {
                active.remove(phase);
                phaseState.phaseSuppression.add(phase);
                primary[static_cast<std::size_t>(sIndex)] = 0.0;
                removed = true;
            }
        }

        if (active.empty())
        {
            // If Newton overshot all three saturations simultaneously, choose
            // the phase whose pre-reset saturation was least negative.
            const auto it = std::max_element(newtonSaturation.begin(), newtonSaturation.end());
            const auto survivor = phases_[static_cast<std::size_t>(std::distance(newtonSaturation.begin(), it))];
            active.add(survivor);
            phaseState.phaseSuppression.remove(survivor);
            primary[static_cast<std::size_t>(three_phase_detail::saturationIndex(
                survivor,
                Indices::Primary::liquidSaturation,
                Indices::Primary::vaporSaturation,
                Indices::Primary::waterSaturation))] = 1.0;
            removed = true;
        }

        updateResult.phaseRemoved = removed;
        phaseState.phasePresence = active;
        phaseState.overallComposition = transitionOverall;

        const auto compositions = normalizedPhaseCompositions_(primary);

        // A disappearing phase first performs a *restricted* equilibrium
        // transform using only the surviving phase identities.  Stability must
        // be evaluated around that reduced equilibrium state; otherwise an
        // unrestricted three-phase flash could immediately recreate a phase
        // that the stability test says is legitimately absent.
        if (removed)
        {
            const auto reduced = flash_.flashRestricted(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                transitionOverall,
                active,
                compositions);

            if (reduced.converged)
            {
                bool unstableMissingPhase = false;
                if (reduced.presence.count() < 3)
                {
                    const auto stability = flash_.stabilityTest(
                        primary[Indices::Primary::pressure],
                        fluid_.temperature,
                        transitionOverall,
                        reduced.presence,
                        reduced.composition);
                    if (!stability.valid)
                    {
                        updateResult.stabilityInvalid = true;
                        unstableMissingPhase = true;
                    }
                    else
                    {
                        unstableMissingPhase =
                            hasStrongMissingPhaseInstability_(
                                stability, reduced.presence,
                                phaseState.phaseSuppression);
                        updateResult.missingPhaseUnstable = unstableMissingPhase;
                    }
                }

                if (!unstableMissingPhase)
                {
                    if (assignCanonicalSingleNonaqueousReducedSet_(
                            primary,
                            phaseState,
                            transitionOverall,
                            reduced.presence,
                            reduced.composition))
                    {
                        updateResult.status = PhaseUpdateStatus::StableReducedSet;
                        updateResult.phaseRemoved = true;
                        return updateResult;
                    }
                    assignFlashResult(primary, phaseState, reduced);
                    updateResult.status = PhaseUpdateStatus::StableReducedSet;
                    updateResult.phaseRemoved = true;
                    return updateResult;
                }
            }

            // A missing phase is unstable (or the restricted transform failed):
            // release the active-set restriction and let the global P-T-z flash
            // select the thermodynamically admissible one/two/three-phase set.
            const auto reflashed = flash_.flash(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                transitionOverall);
            if (reflashed.converged)
            {
                assignReappearingFlashResult_(
                    primary, phaseState, active, reflashed);
                updateResult.phaseRemoved = true;
                return updateResult;
            }

            updateResult.phaseRemoved = true;
            updateResult.restrictedFlashFailed = !reduced.converged;
            updateResult.unrestrictedFlashFailed = true;
            primary = originalPrimary;
            phaseState = originalPhaseState;
            updateResult.status = PhaseUpdateStatus::RecoverableThermodynamicFailure;
            return updateResult;
        }
        else if (active.count() < 3)
        {
            const auto stability = flash_.stabilityTest(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                transitionOverall,
                active,
                compositions);
            const bool stabilityInvalid = !stability.valid;
            const bool missingPhaseUnstable =
                !stabilityInvalid &&
                hasStrongMissingPhaseInstability_(
                    stability, active, phaseState.phaseSuppression);
            if (stabilityInvalid || missingPhaseUnstable)
            {
                updateResult.stabilityInvalid = stabilityInvalid;
                updateResult.missingPhaseUnstable = missingPhaseUnstable;
                const auto reflashed = flash_.flash(
                    primary[Indices::Primary::pressure],
                    fluid_.temperature,
                    transitionOverall);
                if (reflashed.converged)
                {
                    assignReappearingFlashResult_(
                        primary, phaseState, active, reflashed);
                    return updateResult;
                }

                updateResult.unrestrictedFlashFailed = true;
                primary = originalPrimary;
                phaseState = originalPhaseState;
                updateResult.status = PhaseUpdateStatus::RecoverableThermodynamicFailure;
                return updateResult;
            }
        }

        if (active.count() < 3 &&
            assignCanonicalSingleNonaqueousReducedSet_(
                primary,
                phaseState,
                transitionOverall,
                active,
                compositions))
        {
            updateResult.status = PhaseUpdateStatus::StableReducedSet;
            return updateResult;
        }

        // Inactive phases stay numerically finite but carry zero saturation.
        for (CompositionalPhase phase : phases_)
        {
            if (phaseState.phasePresence.contains(phase))
                continue;
            const int sIndex = three_phase_detail::saturationIndex(
                phase,
                Indices::Primary::liquidSaturation,
                Indices::Primary::vaporSaturation,
                Indices::Primary::waterSaturation);
            primary[static_cast<std::size_t>(sIndex)] = 0.0;
            three_phase_detail::writeCompositionToPrimary<Indices>(
                primary,
                transitionOverall,
                phase);
        }
        updateSecondary(primary, phaseState);
        if (phaseState.phasePresence.count() < 3)
            updateResult.status = PhaseUpdateStatus::StableReducedSet;
        return updateResult;
    }

    /**
     * @brief 相态切换前不预先裁剪三相饱和度，保留 active-set 相消失判据所需的信息。
     *
     * 饱和度保持原值直到 `updatePhaseState()`；该函数既可识别负饱和度，
     * 也可识别低于 `phaseBoundaryProbeSaturation` 的边界相。这里若提前
     * clamp，会掩盖 Newton 正在逼近相边界的事实。
     */
    void sanitizePrimaryBeforeFlash(
        PrimaryArray &primary,
        bool useVariableBounds) const
    {
        if (useVariableBounds)
            return;
        // Only protect grossly invalid independent mole fractions; saturation
        // signs are intentionally left untouched.
        for (const auto *indices : {
                 &Indices::Primary::liquidComposition,
                 &Indices::Primary::vaporComposition,
                 &Indices::Primary::waterComposition})
        {
            double independentSum = 0.0;
            for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
            {
                const int index = (*indices)[static_cast<std::size_t>(i)];
                double &value = primary[static_cast<std::size_t>(index)];
                value = std::clamp(value, flash_.options().compositionFloor,
                                   1.0 - flash_.options().compositionFloor);
                independentSum += value;
            }

            // The dependent last component is 1-sum(x_1..x_{N-1}).  Preserve
            // a strictly positive last component instead of allowing individual
            // clamps to create sum(x_independent)>1.
            const double maximumIndependentSum =
                1.0 - flash_.options().compositionFloor;
            if (independentSum > maximumIndependentSum)
            {
                const double scale = maximumIndependentSum / independentSum;
                for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
                {
                    const int index = (*indices)[static_cast<std::size_t>(i)];
                    primary[static_cast<std::size_t>(index)] *= scale;
                }
            }
        }
    }

private:
    /**
     * @brief 保持 SW Newton 未知量采用公开的富油相与富气相顺序。
     *
     * In a one-root nonaqueous region, exchanging the complete Oil and Gas
     * records leaves the SW fugacity equations unchanged.  Newton can therefore
     * converge to either permutation in different cells.  The transport model
     * is not phase symmetric, so canonicalize the accepted state with Wilson
     * volatility ordering before secondary properties and mobilities are used.
     */
    void canonicalizeSwNonAqueousRoles_(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState) const
    {
        if (!fluid_.eos.usesSoreideWhitson() ||
            !phaseState.phasePresence.contains(CompositionalPhase::Oil) ||
            !phaseState.phasePresence.contains(CompositionalPhase::Gas))
        {
            return;
        }

        const auto composition = normalizedPhaseCompositions_(primary);
        const double pressure = primary[Indices::Primary::pressure];
        double orientation = 0.0;
        double magnitude = 0.0;
        for (int component = 0; component < N; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            const double tc = fluid_.eos.mixture().criticalTemperature(component);
            const double pc = fluid_.eos.mixture().criticalPressure(component);
            const double omega = fluid_.eos.mixture().acentricFactor(component);
            const double lnK = std::log(pc / pressure) +
                5.373 * (1.0 + omega) * (1.0 - tc / fluid_.temperature);
            const double term =
                (composition[1][c] - composition[0][c]) * lnK;
            orientation += term;
            magnitude += std::abs(term);
        }

        constexpr double relativeOrientationTolerance = 1.0e-10;
        if (orientation >=
            -relativeOrientationTolerance * std::max(1.0, magnitude))
        {
            return;
        }

        std::swap(primary[Indices::Primary::liquidSaturation],
                  primary[Indices::Primary::vaporSaturation]);
        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            std::swap(
                primary[static_cast<std::size_t>(
                    Indices::Primary::liquidComposition[c])],
                primary[static_cast<std::size_t>(
                    Indices::Primary::vaporComposition[c])]);
        }
        std::swap(phaseState.phaseMoleFraction[0],
                  phaseState.phaseMoleFraction[1]);
        std::swap(phaseState.phaseCompressibility[0],
                  phaseState.phaseCompressibility[1]);

        const auto canonical = normalizedPhaseCompositions_(primary);
        const double floor = flash_.options().compositionFloor;
        for (int component = 0; component < N; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            phaseState.vaporOilEquilibriumRatio[c] =
                canonical[1][c] / std::max(canonical[0][c], floor);
            if (phaseState.phasePresence.contains(CompositionalPhase::Water))
            {
                phaseState.waterOilEquilibriumRatio[c] =
                    canonical[2][c] / std::max(canonical[0][c], floor);
            }
        }
    }

    [[nodiscard]] CompositionalPhase expectedSingleNonAqueousRole_(
        double pressure,
        const Composition &composition) const
    {
        double volatility = 0.0;
        double magnitude = 0.0;
        for (int component = 0; component < N; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            const double tc = fluid_.eos.mixture().criticalTemperature(component);
            const double pc = fluid_.eos.mixture().criticalPressure(component);
            const double omega = fluid_.eos.mixture().acentricFactor(component);
            const double lnK = std::log(pc / pressure) +
                5.373 * (1.0 + omega) * (1.0 - tc / fluid_.temperature);
            const double term = composition[c] * lnK;
            volatility += term;
            magnitude += std::abs(term);
        }

        constexpr double relativeVolatilityTolerance = 1.0e-10;
        return volatility >=
                -relativeVolatilityTolerance * std::max(1.0, magnitude)
            ? CompositionalPhase::Gas
            : CompositionalPhase::Oil;
    }

    [[nodiscard]] bool assignCanonicalSingleNonaqueousReducedSet_(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState,
        const Composition &overallComposition,
        PhasePresence active,
        std::array<Composition, 3> phaseCompositionSeed) const
    {
        if (fluid_.eos.usesSoreideWhitson())
            return false;

        const bool hasOil = active.contains(CompositionalPhase::Oil);
        const bool hasGas = active.contains(CompositionalPhase::Gas);
        if (hasOil == hasGas)
            return false;

        const CompositionalPhase activeNonaqueous =
            hasOil ? CompositionalPhase::Oil : CompositionalPhase::Gas;
        const std::size_t activeIndex =
            static_cast<std::size_t>(phaseIndex(activeNonaqueous));
        const CompositionalPhase expected =
            expectedSingleNonAqueousRole_(
                primary[Indices::Primary::pressure],
                phaseCompositionSeed[activeIndex]);
        if (expected == activeNonaqueous)
            return false;

        PhasePresence corrected = active;
        corrected.remove(activeNonaqueous);
        corrected.add(expected);
        const std::size_t expectedIndex =
            static_cast<std::size_t>(phaseIndex(expected));
        phaseCompositionSeed[expectedIndex] =
            phaseCompositionSeed[activeIndex];

        const auto canonical = flash_.flashRestricted(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            overallComposition,
            corrected,
            phaseCompositionSeed);
        if (!canonical.converged || canonical.presence.bits() != corrected.bits())
            return false;

        assignFlashResult(primary, phaseState, canonical);
        return true;
    }

    /**
     * @brief Assign a global flash while preserving a trace-phase appearance hold.
     *
     * A certified missing-phase instability may yield an equilibrium phase with
     * a physically positive but extremely small saturation.  Such a phase must
     * remain active long enough for the Newton equations to continue on that
     * branch; deleting it again solely because S is still inside the active-set
     * probe band creates deterministic O/W <-> O chatter.  The hold is carried
     * in phaseSuppression until S grows above the probe band.  A zero/negative
     * saturation is never protected and can still disappear normally.
     */
    void assignReappearingFlashResult_(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState,
        PhasePresence previousActive,
        const FlashResult &flash) const
    {
        assignFlashResult(primary, phaseState, flash);
        for (CompositionalPhase phase : phases_)
        {
            if (previousActive.contains(phase) ||
                !flash.presence.contains(phase))
            {
                continue;
            }

            const int sIndex = three_phase_detail::saturationIndex(
                phase,
                Indices::Primary::liquidSaturation,
                Indices::Primary::vaporSaturation,
                Indices::Primary::waterSaturation);
            const double saturation =
                primary[static_cast<std::size_t>(sIndex)];
            if (saturation > 0.0 &&
                saturation <= NaturalNumerics::phaseBoundaryProbeSaturation)
            {
                phaseState.phaseSuppression.add(phase);
            }
        }
    }

    /** @brief 缺失相只有明显越过稳定性边界时才重新生成，避免边界 active-set 抖动。 */
    [[nodiscard]] static bool hasStrongMissingPhaseInstability_(
        const StabilityResult &stability,
        PhasePresence present,
        PhasePresence suppression)
    {
        // invalid stability is handled explicitly by the caller so diagnostics can
        // distinguish "unknown" from a certified unstable missing phase.
        if (!stability.valid)
            return false;

        for (CompositionalPhase phase : phases_)
        {
            if (present.contains(phase))
                continue;
            const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
            const double stabilityMargin = suppression.contains(phase)
                ? NaturalNumerics::phaseHysteresisReappearanceMargin
                : NaturalNumerics::phaseAppearanceStabilityMargin;
            if (stability.missingPhaseUnstable[p] &&
                stability.trialSum[p] >
                    1.0 + stabilityMargin)
                return true;
        }
        return false;
    }

    [[nodiscard]] static ThreePhaseFlashOptions flashOptions_(const FluidSystem<Indices> &fluid)
    {
        ThreePhaseFlashOptions options;
        options.waterComponent = fluid.fullyCompositionalWaterComponent;
        return options;
    }

    /** @brief 从 Natural 主变量读取并统一归一化 O/G/W 三相组成。 */
    [[nodiscard]] std::array<Composition, 3> normalizedPhaseCompositions_(
        const PrimaryArray &primary) const
    {
        std::array<Composition, 3> composition{};
        const double floor = flash_.options().compositionFloor;
        for (CompositionalPhase phase : phases_)
        {
            const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
            composition[p] =
                three_phase_detail::compositionFromPrimary<Indices>(primary, phase);
            three_phase_detail::normalize(composition[p], floor);
        }
        return composition;
    }

    [[nodiscard]] Composition overallFromFlash_(const FlashResult &flash) const
    {
        Composition z{};
        for (int p = 0; p < 3; ++p)
        {
            for (int i = 0; i < N; ++i)
            {
                z[static_cast<std::size_t>(i)] +=
                    flash.phaseMoleFraction[static_cast<std::size_t>(p)] *
                    flash.composition[static_cast<std::size_t>(p)][static_cast<std::size_t>(i)];
            }
        }
        three_phase_detail::normalize(z, flash_.options().compositionFloor);
        return z;
    }

    [[nodiscard]] Composition transitionOverallComposition_(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        const auto composition = normalizedPhaseCompositions_(primary);
        Composition z{};
        double betaSum = 0.0;

        for (CompositionalPhase phase : phases_)
        {
            const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
            if (!phaseState.phasePresence.contains(phase))
                continue;
            const double beta = std::max(phaseState.phaseMoleFraction[p], 0.0);
            betaSum += beta;
            for (int i = 0; i < N; ++i)
                z[static_cast<std::size_t>(i)] += beta * composition[p][static_cast<std::size_t>(i)];
        }

        if (!(betaSum > NaturalNumerics::minimumNormalizationDenominator))
            return phaseState.overallComposition;

        for (double &value : z) value /= betaSum;
        three_phase_detail::normalize(z, flash_.options().compositionFloor);
        return z;
    }

    [[nodiscard]] static CompositionalPhase firstActive_(PhasePresence active)
    {
        for (CompositionalPhase phase : phases_)
            if (active.contains(phase)) return phase;
        throw std::logic_error("No active phase.");
    }

    inline static constexpr std::array<CompositionalPhase, 3> phases_{
        CompositionalPhase::Oil,
        CompositionalPhase::Gas,
        CompositionalPhase::Water};

    const FluidSystem<Indices> &fluid_;
    Flash flash_;
};

} // namespace MPMC
