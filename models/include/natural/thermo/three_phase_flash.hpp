/**
 * @file three_phase_flash.hpp
 * @brief O/G/W 三相 flash、稳定性测试及受证书约束的失败恢复算法。
 */
#pragma once

#include <common/math.hpp>
#include <natural/numerics.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/three_phase_flash_types.hpp>
#include <natural/thermo/phase_behavior.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace MPMC
{

/**
 * @brief 与 EOS 类型解耦的 O/G/W 三相 flash 和多相切平面稳定性测试。
 *
 * 公共槽位约定：0=油富集液相，1=气相，2=水富集液相；PR/SW/CPA 负责各自的根选择和逸度。
 * 所有相使用同一组分集，因此 H2O 可进入油/气相，烃和 CO2 也可进入水富集相。
 *
 * 数学：以油相为参考，`K_gi=y_i/x_i`、`K_wi=a_i/x_i`，
 * `beta_o=1-beta_g-beta_w`，`D_i=beta_o+beta_g K_gi+beta_w K_wi`，`x_i=z_i/D_i`。
 * 该 `D_i` 写法与传统式完全等价，但在 `beta_o` 很小或 K 值跨越多数量级时可减少消减误差。
 * 广义 Rachford–Rice 方程分别令 `sum z_i(K_gi-1)/D_i=0` 和
 * `sum z_i(K_wi-1)/D_i=0`；逸度相等通过 `K_g=phi_o/phi_g`、
 * `K_w=phi_o/phi_w` 更新。相摩尔分率最终按 EOS 摩尔体积换算为孔隙饱和度。
 */
template <class Indices>
class CubicThreePhaseFlash final
{
public:
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;
    using Result = ThreePhaseFlashResult<Indices>;
    using StabilityResult = ThreePhaseStabilityResult<Indices>;
    using Eos = CubicEquationOfState<Indices>;

    CubicThreePhaseFlash(
        const Eos &eos,
        ThreePhaseFlashOptions options = {})
        : eos_(eos), options_(options)
    {
        if (options_.waterComponent >= N)
            throw std::invalid_argument("Three-phase flash water-component index is out of range.");
        if (options_.maximumIterations <= 0 ||
            options_.maximumRachfordRiceIterations <= 0 ||
            options_.maximumStabilityIterations <= 0 ||
            options_.stabilitySsiIterationsBeforeNewton < 0 ||
            options_.maximumStabilityNewtonIterations < 0 ||
            options_.twoPhaseSsiIterationsBeforeNewton < 0 ||
            options_.maximumTwoPhaseNewtonIterations < 0 ||
            options_.threePhaseSsiIterationsBeforeNewton < 0 ||
            options_.maximumThreePhaseNewtonIterations < 0)
        {
            throw std::invalid_argument(
                "Three-phase flash iteration limits must be positive/non-negative.");
        }
        if (!(options_.fugacityTolerance > 0.0) ||
            !(options_.rachfordRiceTolerance > 0.0) ||
            !(options_.stabilityTolerance > 0.0) ||
            !(options_.phaseFractionTolerance >= 0.0) ||
            !(options_.coincidentPhaseCompositionTolerance >= 0.0) ||
            !(options_.compositionFloor > 0.0) ||
            !(options_.logKStepLimit > 0.0) ||
            !(options_.numericalJacobianStep > 0.0))
        {
            throw std::invalid_argument(
                "Three-phase flash tolerances and numerical steps must be valid.");
        }
    }

    [[nodiscard]] const ThreePhaseFlashOptions &options() const noexcept { return options_; }

    /** @brief 由压力 [Pa]、温度 [K] 和总摩尔分数 `z` 计算稳定相态及各相组成。 */
    [[nodiscard]] Result flash(
        double pressure,
        double temperature,
        Composition z) const
    {
        validatePT_(pressure, temperature);
        normalize_(z);

        Composition kg = wilsonGasOilK_(pressure, temperature);
        Composition kw = initialWaterOilK_(z);

        // 数值：相数由“当前活动相集平衡 + 缺失相稳定性”决定，不能把第一个收敛的三相代数根
        // 当成物理解。先从常见 O+G（不可行则单相）出发，仅在 TPD 表明需要时才生成缺失相。
        Result seed = flashTwo_(pressure, temperature, z,
                                CompositionalPhase::Oil,
                                CompositionalPhase::Gas,
                                kg);
        PhasePresence active = presenceFromFractions_(seed.phaseMoleFraction);
        if (!seed.converged || active.empty())
        {
            active = chooseSinglePhase_(pressure, temperature, z);
            seed = singlePhaseResult_(pressure, temperature, z, active);
        }

        Result selected = equilibrateActiveSet_(
            pressure, temperature, z, active, seed);
        if (selected.converged)
        {
            canonicalizeFinal_(pressure, temperature, selected);
            if (selected.converged &&
                validatePhysicalFlashResult_(pressure, temperature, z, selected))
            {
                Result preferred = preferStableSwReducedOverSingle_(
                    pressure, temperature, z, selected);
                if (preferred.converged)
                    return preferred;
            }
            selected.converged = false;
        }

        // 数值：直接三相求解只作为 active-set 路径失败后的恢复手段，不承担相数判定职责。
        Result full = flashThree_(pressure, temperature, z, kg, kw);
        if (full.converged)
        {
            canonicalize_(pressure, temperature, full);
            canonicalizeFinal_(pressure, temperature, full);
            if (full.converged &&
                validatePhysicalFlashResult_(pressure, temperature, z, full))
            {
                Result preferred = preferStableSwReducedOverSingle_(
                    pressure, temperature, z, full);
                if (preferred.converged)
                    return preferred;
            }
            full.converged = false;
        }

        // 数值：active-set 与直接三相路径都失败时，先对所有 EOS 枚举稳定的一/两相
        // 约化子集。该搜索只接受 stability.valid && stability.stable 的候选并按 Gibbs
        // 最小值选择，因此它是通用热力学 fallback，而不是 SW 专属数值补丁。
        Result reduced = stableReducedSetFallback_(
            pressure, temperature, z);
        if (reduced.converged)
        {
            canonicalizeFinal_(pressure, temperature, reduced);
            if (reduced.converged &&
                validatePhysicalFlashResult_(pressure, temperature, z, reduced))
                return reduced;
        }

        // Exact-mass allocation is a fail-only basin search. SW keeps
        // its historical use. The Jia external CPA audit may opt in after
        // ordinary active-set/direct/reduced paths have all failed; active
        // phase root identities remain fixed by phaseResult_ and no EOS
        // parameter/stability criterion is changed. Production CPA keeps the
        // option false and therefore does not enter this recovery layer.
        const bool useAllocationRecovery =
            eos_.usesSoreideWhitson() ||
            (eos_.usesCubicPlusAssociation() &&
             options_.cpaSelectGibbsMinimumRoot);
        if (useAllocationRecovery)
        {
            Result allocated = flashThreeAllocationFallback_(
                pressure, temperature, z);
            if (allocated.converged)
            {
                canonicalizeFinal_(pressure, temperature, allocated);
                if (allocated.converged &&
                    validatePhysicalFlashResult_(pressure, temperature, z, allocated))
                {
                    Result preferred = preferStableSwReducedOverSingle_(
                        pressure, temperature, z, allocated);
                    if (preferred.converged)
                        return preferred;
                }
            }
        }

        // 数值：SW 的内部热力学角色与公开 O/G/W 槽位解耦。若默认角色分配不能得到自洽平衡，
        // 枚举其余固定角色排列；一次非线性求解期间 Aqueous/Liquid/Vapor 角色保持不变，
        // 因而不会仅因输出槽位重命名而切换水相 BIP 或液/气根。
        if (eos_.usesSoreideWhitson() && enableSwRolePermutationSearch_)
        {
            Result bestPermuted;
            double bestG = std::numeric_limits<double>::infinity();
            for (const auto &roles : swRolePermutations_)
            {
                if (roles == thermodynamicRoles_)
                    continue;
                CubicThreePhaseFlash candidate(
                    eos_, options_, roles, false, false);
                Result permuted = candidate.flash(pressure, temperature, z);
                if (!permuted.converged)
                    continue;
                const double g = candidate.dimensionlessGibbs_(
                    pressure, temperature, permuted);
                if (std::isfinite(g) && g < bestG)
                {
                    bestG = g;
                    bestPermuted = permuted;
                }
            }
            if (bestPermuted.converged)
            {
                Result preferred = preferStableSwReducedOverSingle_(
                    pressure, temperature, z, bestPermuted);
                if (preferred.converged)
                    return preferred;
            }
        }

        // 数值：压力延拓严格放在全部固定角色路径之后，保持已有收敛状态的结果不变。
        // PR/SW/CPA 都只在全部常规/恢复分支失败时进入同伦，因此不会改变正常收敛路径。
        if (enablePressureContinuation_)
        {
            Result continued = pressureContinuationFallback_(
                pressure, temperature, z);
            if (continued.converged)
            {
                canonicalizeFinal_(pressure, temperature, continued);
                if (continued.converged &&
                    validatePhysicalFlashResult_(pressure, temperature, z, continued))
                {
                    Result preferred = preferStableSwReducedOverSingle_(
                        pressure, temperature, z, continued);
                    if (preferred.converged)
                        return preferred;
                }
                continued.converged = false;
            }
        }
        return selected;
    }

    /**
     * @brief 在禁止 `allowed` 之外相出现的条件下重新平衡 P–T–z。
     *
     * 状态：某相因负饱和度被移除后使用该约化变换；若剩余相分率也趋零，可继续约化。
     * 本函数不会生成新相；新相只能由 `stabilityTest()` 判定后再进入 unrestricted `flash()`。
     */
    [[nodiscard]] Result flashRestricted(
        double pressure,
        double temperature,
        Composition z,
        PhasePresence allowed) const
    {
        normalize_(z);
        Result result = flashRestrictedImpl_(
            pressure, temperature, z, allowed, nullptr);
        if (!result.converged && allowed.count() == 2 &&
            eos_.usesSoreideWhitson())
        {
            const auto pair = activePair_(allowed);
            finalizePhysicallyConvergedTwoPhase_(
                pressure, temperature, z, pair.first, pair.second, result);
        }
        if (result.converged)
        {
            canonicalizeFinal_(pressure, temperature, result, true);
            if (result.converged &&
                !validatePhysicalFlashResult_(pressure, temperature, z, result))
                result.converged = false;
        }
        return result;
    }

    /**
     * @brief 使用当前各相组成作为初值的 restricted flash。
     *
     * 数值：气相刚消失转入 O+W 时，切换前组成通常比通用水/油 K 初值更接近正确液液平衡。
     */
    [[nodiscard]] Result flashRestricted(
        double pressure,
        double temperature,
        Composition z,
        PhasePresence allowed,
        const std::array<Composition, 3> &phaseCompositionSeed) const
    {
        normalize_(z);
        Result result = flashRestrictedImpl_(
            pressure, temperature, z, allowed,
            &phaseCompositionSeed);
        if (!result.converged && allowed.count() == 2 &&
            eos_.usesSoreideWhitson())
        {
            const auto pair = activePair_(allowed);
            finalizePhysicallyConvergedTwoPhase_(
                pressure, temperature, z, pair.first, pair.second, result);
        }
        if (result.converged)
        {
            canonicalizeFinal_(pressure, temperature, result, true);
            if (result.converged &&
                !validatePhysicalFlashResult_(pressure, temperature, z, result))
                result.converged = false;
        }
        return result;
    }

    /**
     * @brief 对当前一/两相平衡逐一执行缺失相切平面稳定性测试。
     *
     * 物理：精确多相平衡时各活动相具有相同组分化学势，因此任选一个活动相作为参考。
     * 数值：每个缺失相尝试多组相特定初值，并保留最不稳定的驻点，避免单一初值漏检新相。
     */
    [[nodiscard]] StabilityResult stabilityTest(
        double pressure,
        double temperature,
        const Composition &zInput,
        PhasePresence active,
        const std::array<Composition, 3> &activeCompositions) const
    {
        validatePT_(pressure, temperature);
        if (active.empty())
            throw std::invalid_argument("Stability test requires at least one active phase.");

        Composition z = zInput;
        normalize_(z);
        StabilityResult result;
        result.testedPresence = active;

        const CompositionalPhase referencePhase = firstActive_(active);
        const int ref = phaseIndex(referencePhase);
        std::array<Composition, 3> stabilityCompositions = activeCompositions;

        // A single active phase must carry the overall composition exactly.
        // For SW this also prevents a stale inactive-slot composition from
        // becoming the reference chemical potential after an O/G role switch.
        if (eos_.usesSoreideWhitson() && active.count() == 1)
            stabilityCompositions[static_cast<std::size_t>(ref)] = z;

        for (CompositionalPhase phase : phases_)
        {
            if (active.contains(phase))
                normalize_(stabilityCompositions[static_cast<std::size_t>(phaseIndex(phase))]);
        }

        // SW assigns different root/BIP semantics to nonaqueous and aqueous
        // roles.  Natural calls stabilityTest() while Newton is still moving,
        // so the supplied active-phase compositions are not necessarily at
        // interphase fugacity equilibrium yet.  Building TPD directly from
        // such an iterate can turn the existing aqueous basin into a false
        // "missing Gas" certificate.  Only when the active records fail a
        // cheap fugacity-closure check, first certify the current active set
        // with a restricted flash; already-equilibrated states pay no extra
        // flash cost.
        if (eos_.usesSoreideWhitson() && active.count() > 1)
        {
            Composition rawReference =
                stabilityCompositions[static_cast<std::size_t>(ref)];
            const auto rawReferenceThermo = phaseResult_(
                referencePhase, pressure, temperature, rawReference);
            double activeFugacityMismatch = 0.0;
            for (CompositionalPhase phase : phases_)
            {
                if (!active.contains(phase) || phase == referencePhase)
                    continue;
                const std::size_t slot =
                    static_cast<std::size_t>(phaseIndex(phase));
                const auto thermo = phaseResult_(
                    phase, pressure, temperature, stabilityCompositions[slot]);
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(i);
                    const double xRef = rawReference[c];
                    const double x = stabilityCompositions[slot][c];
                    if (xRef <= 10.0 * options_.compositionFloor ||
                        x <= 10.0 * options_.compositionFloor)
                        continue;
                    const double refActivity = std::max(
                        xRef * rawReferenceThermo.fugacityCoefficient[c],
                        1.0e-300);
                    const double activity = std::max(
                        x * thermo.fugacityCoefficient[c], 1.0e-300);
                    activeFugacityMismatch = std::max(
                        activeFugacityMismatch,
                        std::abs(std::log(refActivity / activity)));
                }
            }

            const double certificationTolerance = std::max(
                10.0 * options_.fugacityTolerance,
                10.0 * options_.stabilityTolerance);
            if (activeFugacityMismatch > certificationTolerance)
            {
                Result certified = flashRestrictedImpl_(
                    pressure, temperature, z, active, &stabilityCompositions);
                if (!certified.converged ||
                    certified.presence.bits() != active.bits())
                {
                    // The supplied active set is itself not a certified
                    // equilibrium state.  Do not manufacture a missing-phase
                    // TPD decision from it; the caller may invoke its existing
                    // full-flash recovery path.
                    result.valid = false;
                    result.stable = false;
                    return result;
                }
                stabilityCompositions = certified.composition;
            }
        }

        Composition referenceComposition =
            stabilityCompositions[static_cast<std::size_t>(ref)];
        normalize_(referenceComposition);
        const auto referenceThermo = phaseResult_(
            referencePhase, pressure, temperature, referenceComposition);

        std::array<double, N> logReference{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double xi = std::max(referenceComposition[c], options_.compositionFloor);
            const double phi = std::max(referenceThermo.fugacityCoefficient[c], 1.0e-300);
            logReference[c] = std::log(xi) + std::log(phi);
        }

        for (CompositionalPhase candidate : phases_)
        {
            if (active.contains(candidate))
                continue;

            const int p = phaseIndex(candidate);
            Trial best;
            best.sum = -std::numeric_limits<double>::infinity();

            const auto seeds = stabilitySeeds_(
                candidate, pressure, temperature, z, referenceComposition);
            for (std::size_t seedIndex = 0; seedIndex < seeds.size; ++seedIndex)
            {
                const auto &seed = seeds.values[seedIndex];
                const Trial trial = stabilityTrial_(
                    candidate, pressure, temperature, logReference,
                    referenceComposition, seed);
                if (trial.valid && (!best.valid || trial.sum > best.sum))
                    best = trial;
            }

            if (!best.valid)
            {
                result.valid = false;
                result.stable = false; // 数值：试探失败时按“不稳定”保守处理，避免错误抑制新相出现。
                continue;
            }

            result.trialSum[static_cast<std::size_t>(p)] = best.sum;
            result.incipientComposition[static_cast<std::size_t>(p)] = best.composition;

            // A TPD solver is allowed to converge to a stationary point that
            // is simply one of the phases already present.  That point is
            // trivial even when it is not the first/reference phase.  Without
            // this check a Water-like active phase can be rediscovered through
            // the SW Gas/Oil root and be misreported as a new phase.
            bool coincidesWithActivePhase = false;
            const double coincidenceTolerance = std::max(
                10.0 * options_.coincidentPhaseCompositionTolerance, 1.0e-7);
            for (CompositionalPhase phase : phases_)
            {
                if (!active.contains(phase))
                    continue;
                const std::size_t slot =
                    static_cast<std::size_t>(phaseIndex(phase));
                double maxDifference = 0.0;
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(i);
                    maxDifference = std::max(
                        maxDifference,
                        std::abs(best.composition[c] -
                                 stabilityCompositions[slot][c]));
                }
                if (maxDifference <= coincidenceTolerance)
                {
                    coincidesWithActivePhase = true;
                    break;
                }
            }

            bool unstable = !best.trivial && !coincidesWithActivePhase &&
                best.sum > 1.0 + options_.stabilityTolerance;
            if (candidate == CompositionalPhase::Water && unstable &&
                !eos_.aqueousVolumeCompositionSupported(best.composition))
            {
                // The public Water role selects the configured aqueous volume
                // and viscosity closures.  A TPD direction outside that
                // explicit composition domain is not an admissible phase in
                // this flow model.  Retain its trial sum/composition for
                // diagnostics, but do not release it into the active set.
                unstable = false;
            }
            result.missingPhaseUnstable[static_cast<std::size_t>(p)] = unstable;
            result.stable = result.stable && !unstable;
        }
        return result;
    }

    struct ThreePhaseRachfordRiceResult
    {
        bool converged{false};
        double betaOil{0.0};
        double betaGas{0.0};
        double betaWater{0.0};
    };

    /**
     * @brief 求解广义三相 Rachford-Rice 方程。
     *
     * This small public diagnostic hook is intentionally EOS-independent.  It
     * is useful for validating the multiphase material-balance kernel against
     * published flash tables without conflating that check with fugacity-model
     * differences (PR vs. SW vs. CPA).
     */
    [[nodiscard]] ThreePhaseRachfordRiceResult solveThreePhaseRachfordRice(
        Composition z,
        const Composition &gasOilK,
        const Composition &waterOilK,
        double initialBetaGas = 0.2,
        double initialBetaWater = 0.2) const
    {
        normalize_(z);
        const auto rr = solveRr3Robust_(
            z, gasOilK, waterOilK, initialBetaGas, initialBetaWater);
        if (!rr.converged)
            return {};
        return {
            true,
            1.0 - rr.betaGas - rr.betaWater,
            rr.betaGas,
            rr.betaWater};
    }

    /** @brief 将相摩尔分率与 EOS 压缩因子 Z 转换为饱和度。 */
    void updateSaturationsFromMoles(
        double pressure,
        double temperature,
        Result &result) const
    {
        double volumeSum = 0.0;
        std::array<double, 3> volume{};
        for (int p = 0; p < 3; ++p)
        {
            if (!result.presence.contains(phases_[static_cast<std::size_t>(p)]))
            {
                result.phaseMoleFraction[static_cast<std::size_t>(p)] = 0.0;
                result.saturation[static_cast<std::size_t>(p)] = 0.0;
                result.molarDensity[static_cast<std::size_t>(p)] = 0.0;
                continue;
            }
            const double Z = result.compressibility[static_cast<std::size_t>(p)];
            if (!(Z > 0.0))
                throw std::runtime_error("Three-phase flash produced non-positive compressibility.");
            const auto role = thermodynamicRoles_[static_cast<std::size_t>(p)];
            if (role == CompositionalPhase::Water &&
                !eos_.aqueousVolumeCompositionSupported(
                    result.composition[static_cast<std::size_t>(p)]))
            {
                // A flash iterate can satisfy algebraic fugacity equations yet
                // lie outside an explicitly configured physical-property
                // closure.  Reject that basin so active-set/role-permutation
                // recovery may continue; it is not a valid converged state.
                result.converged = false;
                result.saturation.fill(0.0);
                result.molarDensity.fill(0.0);
                return;
            }
            const double rhoM = eos_.molarDensity(
                pressure,
                temperature,
                result.composition[static_cast<std::size_t>(p)],
                Z,
                role);
            result.molarDensity[static_cast<std::size_t>(p)] = rhoM;
            volume[static_cast<std::size_t>(p)] =
                std::max(result.phaseMoleFraction[static_cast<std::size_t>(p)], 0.0) / rhoM;
            volumeSum += volume[static_cast<std::size_t>(p)];
        }
        if (!(volumeSum > 0.0))
            throw std::runtime_error("Three-phase flash has zero total phase volume.");
        for (int p = 0; p < 3; ++p)
            result.saturation[static_cast<std::size_t>(p)] = volume[static_cast<std::size_t>(p)] / volumeSum;
    }

private:
#include <natural/thermo/detail/three_phase_flash_core.inc>

#include <natural/thermo/detail/three_phase_flash_nonlinear.inc>

#include <natural/thermo/detail/three_phase_flash_restricted.inc>

#include <natural/thermo/detail/three_phase_flash_allocation.inc>

#include <natural/thermo/detail/three_phase_flash_canonicalization.inc>

#include <natural/thermo/detail/three_phase_flash_stability.inc>

    const Eos &eos_;
    ThreePhaseFlashOptions options_;
    std::array<CompositionalPhase, 3> thermodynamicRoles_{
        CompositionalPhase::Oil,
        CompositionalPhase::Gas,
        CompositionalPhase::Water};
    bool enableSwRolePermutationSearch_{true};
    bool enablePressureContinuation_{true};
    mutable CubicParameterCache cubicParameterCache_{};
};

/**
 * @brief 兼容 v31 以前 flash 类名的别名。
 *
 * The implementation has supported PR, Soreide-Whitson, and CPA backends
 * through CubicEquationOfState for several releases.  CubicThreePhaseFlash is
 * therefore the canonical name; this alias preserves existing user code.
 */
template <class Indices>
using PengRobinsonThreePhaseFlash = CubicThreePhaseFlash<Indices>;

} // namespace MPMC
