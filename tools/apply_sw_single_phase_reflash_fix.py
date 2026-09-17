#!/usr/bin/env python3
from pathlib import Path


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one anchor, found {count}")
    path.write_text(text.replace(old, new, 1))


root = Path(__file__).resolve().parents[1]
equilibrium = root / "models/include/natural/state/three_phase_equilibrium.hpp"
test = root / "test/src/unit/sw_phase_ordering_test.cpp"

old_transition = r'''        else if (active.count() < 3)
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
                    assignFlashResult(primary, phaseState, reflashed);
                    return updateResult;
                }

                updateResult.unrestrictedFlashFailed = true;
                primary = originalPrimary;
                phaseState = originalPhaseState;
                updateResult.status = PhaseUpdateStatus::RecoverableThermodynamicFailure;
                return updateResult;
            }
        }
'''

new_transition = r'''        else if (active.count() < 3)
        {
            const auto stability = flash_.stabilityTest(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                transitionOverall,
                active,
                compositions);
            const bool stabilityInvalid = !stability.valid;
            const PhasePresence stronglyUnstableMissing =
                stabilityInvalid
                    ? emptyPhasePresence_()
                    : stronglyUnstableMissingPhases_(
                          stability, active, phaseState.phaseSuppression);
            const bool missingPhaseUnstable = !stronglyUnstableMissing.empty();
            if (stabilityInvalid || missingPhaseUnstable)
            {
                updateResult.stabilityInvalid = stabilityInvalid;
                updateResult.missingPhaseUnstable = missingPhaseUnstable;

                // SW uses different thermodynamic roles for the aqueous and
                // non-aqueous phase models.  Near a single-phase boundary an
                // unrestricted flash can therefore converge to a different
                // single-role basin than the missing phase that TPD actually
                // certified.  First expand only the certified active set and
                // seed the new phase with the TPD stationary composition.  The
                // expanded state must itself pass a fresh stability certificate.
                if (!stabilityInvalid &&
                    tryAssignSwTriggeredExpandedSet_(
                        primary,
                        phaseState,
                        transitionOverall,
                        active,
                        stability,
                        phaseState.phaseSuppression,
                        compositions))
                {
                    return updateResult;
                }

                const auto reflashed = flash_.flash(
                    primary[Indices::Primary::pressure],
                    fluid_.temperature,
                    transitionOverall);
                if (reflashed.converged &&
                    swUnrestrictedSinglePhaseConsistentWithTrigger_(
                        reflashed, active, stronglyUnstableMissing))
                {
                    assignFlashResult(primary, phaseState, reflashed);
                    return updateResult;
                }

                updateResult.unrestrictedFlashFailed = true;
                primary = originalPrimary;
                phaseState = originalPhaseState;
                updateResult.status = PhaseUpdateStatus::RecoverableThermodynamicFailure;
                return updateResult;
            }
        }
'''
replace_once(equilibrium, old_transition, new_transition)

old_helper = r'''    /** @brief 缺失相只有明显越过稳定性边界时才重新生成，避免边界 active-set 抖动。 */
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
'''

new_helper = r'''    [[nodiscard]] static PhasePresence emptyPhasePresence_()
    {
        PhasePresence empty = PhasePresence::all();
        for (CompositionalPhase phase : phases_)
            empty.remove(phase);
        return empty;
    }

    /** @brief 返回真正越过 appearance/hysteresis margin 的缺失相集合。 */
    [[nodiscard]] static PhasePresence stronglyUnstableMissingPhases_(
        const StabilityResult &stability,
        PhasePresence present,
        PhasePresence suppression)
    {
        PhasePresence unstable = emptyPhasePresence_();
        if (!stability.valid)
            return unstable;

        for (CompositionalPhase phase : phases_)
        {
            if (present.contains(phase))
                continue;
            const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
            const double stabilityMargin = suppression.contains(phase)
                ? NaturalNumerics::phaseHysteresisReappearanceMargin
                : NaturalNumerics::phaseAppearanceStabilityMargin;
            if (stability.missingPhaseUnstable[p] &&
                stability.trialSum[p] > 1.0 + stabilityMargin)
            {
                unstable.add(phase);
            }
        }
        return unstable;
    }

    /** @brief 缺失相只有明显越过稳定性边界时才重新生成，避免边界 active-set 抖动。 */
    [[nodiscard]] static bool hasStrongMissingPhaseInstability_(
        const StabilityResult &stability,
        PhasePresence present,
        PhasePresence suppression)
    {
        return !stronglyUnstableMissingPhases_(
                    stability, present, suppression)
                    .empty();
    }

    /**
     * @brief SW phase appearance 先沿 TPD 证书扩展 active set，再决定是否需要全局 reflash。
     *
     * A single-phase SW state can have several algebraic single-role basins
     * because the Aqueous role changes the H2O BIP model.  Once TPD has certified
     * a particular missing phase, replacing the cell by an unrelated one-phase
     * role is not a continuous phase-appearance operation.  Seed only the
     * certified missing phase(s), solve that restricted set, and require the
     * resulting equilibrium to be stable before committing it.
     */
    [[nodiscard]] bool tryAssignSwTriggeredExpandedSet_(
        PrimaryArray &primary,
        PhaseStateData<Indices> &phaseState,
        const Composition &overallComposition,
        PhasePresence active,
        const StabilityResult &stability,
        PhasePresence suppression,
        std::array<Composition, 3> phaseCompositionSeed) const
    {
        if (!fluid_.eos.usesSoreideWhitson() || !stability.valid)
            return false;

        const PhasePresence triggered = stronglyUnstableMissingPhases_(
            stability, active, suppression);
        if (triggered.empty())
            return false;

        PhasePresence expanded = active;
        for (CompositionalPhase phase : phases_)
        {
            if (!triggered.contains(phase))
                continue;
            expanded.add(phase);
            const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
            phaseCompositionSeed[p] = stability.incipientComposition[p];
        }

        // Three certified phases already mean a genuine unrestricted problem.
        if (expanded.count() >= 3)
            return false;

        const auto candidate = flash_.flashRestricted(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            overallComposition,
            expanded,
            phaseCompositionSeed);
        if (!candidate.converged)
            return false;

        bool carriesTriggeredPhase = false;
        for (CompositionalPhase phase : phases_)
        {
            carriesTriggeredPhase = carriesTriggeredPhase ||
                (triggered.contains(phase) && candidate.presence.contains(phase));
        }
        if (!carriesTriggeredPhase)
            return false;

        if (candidate.presence.count() < 3)
        {
            const auto certificate = flash_.stabilityTest(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                overallComposition,
                candidate.presence,
                candidate.composition);
            if (!certificate.valid ||
                hasStrongMissingPhaseInstability_(
                    certificate, candidate.presence, suppression))
            {
                return false;
            }
        }

        assignFlashResult(primary, phaseState, candidate);
        return true;
    }

    /**
     * @brief 拒绝与当前 TPD 证书无关的 SW 单相 reflash 角色跳变。
     *
     * A valid one-phase replacement must either retain the current active role
     * or be one of the missing roles that actually crossed the appearance
     * margin.  Multiphase results remain unrestricted and invalid-stability
     * recovery keeps the historical behavior.
     */
    [[nodiscard]] bool swUnrestrictedSinglePhaseConsistentWithTrigger_(
        const FlashResult &reflashed,
        PhasePresence active,
        PhasePresence stronglyUnstableMissing) const
    {
        if (!fluid_.eos.usesSoreideWhitson() ||
            reflashed.presence.count() != 1 ||
            stronglyUnstableMissing.empty())
        {
            return true;
        }

        const CompositionalPhase returned = firstActive_(reflashed.presence);
        return active.contains(returned) ||
            stronglyUnstableMissing.contains(returned);
    }
'''
replace_once(equilibrium, old_helper, new_helper)

fake_flash = r'''
class TriggeredSwFlashStub
{
public:
    using Result = MPMC::ThreePhaseFlashResult<Indices>;
    using StabilityResult = MPMC::ThreePhaseStabilityResult<Indices>;

    TriggeredSwFlashStub(const Eos &, MPMC::ThreePhaseFlashOptions options)
        : options_(std::move(options))
    {}

    [[nodiscard]] const MPMC::ThreePhaseFlashOptions &options() const noexcept
    {
        return options_;
    }

    static void resetCounters()
    {
        restrictedCalls = 0;
        unrestrictedCalls = 0;
    }

    [[nodiscard]] Result flash(
        double,
        double,
        Composition z) const
    {
        ++unrestrictedCalls;
        return makeResult_(MPMC::PhasePresence::waterOnly(), z);
    }

    [[nodiscard]] Result flashRestricted(
        double,
        double,
        Composition z,
        MPMC::PhasePresence allowed) const
    {
        ++restrictedCalls;
        return makeResult_(allowed, z);
    }

    [[nodiscard]] Result flashRestricted(
        double,
        double,
        Composition z,
        MPMC::PhasePresence allowed,
        const std::array<Composition, 3> &) const
    {
        ++restrictedCalls;
        return makeResult_(allowed, z);
    }

    [[nodiscard]] StabilityResult stabilityTest(
        double,
        double,
        const Composition &,
        MPMC::PhasePresence active,
        const std::array<Composition, 3> &) const
    {
        StabilityResult result;
        result.testedPresence = active;
        if (active.bits() == MPMC::PhasePresence::oilBit)
        {
            result.stable = false;
            result.missingPhaseUnstable[1] = true;
            result.trialSum[1] = 1.0 + 1.0e-4;
            result.incipientComposition[1] = {0.10, 0.80, 0.10};
        }
        return result;
    }

    inline static int restrictedCalls{0};
    inline static int unrestrictedCalls{0};

private:
    [[nodiscard]] static Result makeResult_(
        MPMC::PhasePresence presence,
        const Composition &z)
    {
        Result result;
        result.converged = true;
        result.presence = presence;
        result.phaseMoleFraction = {0.0, 0.0, 0.0};
        result.saturation = {0.0, 0.0, 0.0};
        result.composition = {z, z, z};
        const double share = 1.0 / static_cast<double>(presence.count());
        for (MPMC::CompositionalPhase phase : {
                 MPMC::CompositionalPhase::Oil,
                 MPMC::CompositionalPhase::Gas,
                 MPMC::CompositionalPhase::Water})
        {
            if (!presence.contains(phase))
                continue;
            const std::size_t p = static_cast<std::size_t>(MPMC::phaseIndex(phase));
            result.phaseMoleFraction[p] = share;
            result.saturation[p] = share;
        }
        result.compressibility = {1.0, 1.0, 1.0};
        result.vaporOilK.fill(1.0);
        result.waterOilK.fill(1.0);
        return result;
    }

    MPMC::ThreePhaseFlashOptions options_;
};
'''
replace_once(
    test,
    "\nvoid writeComposition(\n",
    fake_flash + "\nvoid writeComposition(\n")

new_test = r'''
void checkSwTriggeredReflashPreservesTpdPhaseIdentity()
{
    using StubEquilibrium =
        MPMC::FullyCompositionalThreePhaseEquilibrium<Indices, TriggeredSwFlashStub>;

    auto fluid = makeFluid();
    StubEquilibrium equilibrium(fluid);
    TriggeredSwFlashStub::resetCounters();

    constexpr Composition overall{0.20, 0.35, 0.45};
    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 5.2e6;
    primary[Indices::Primary::liquidSaturation] = 1.0;
    primary[Indices::Primary::vaporSaturation] = 0.0;
    primary[Indices::Primary::waterSaturation] = 0.0;
    writeComposition(primary, Indices::Primary::liquidComposition, overall);
    writeComposition(primary, Indices::Primary::vaporComposition, overall);
    writeComposition(primary, Indices::Primary::waterComposition, overall);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::oilOnly();
    phaseState.phaseMoleFraction = {1.0, 0.0, 0.0};
    phaseState.overallComposition = overall;

    const auto update = equilibrium.updatePhaseState(primary, phaseState);
    const auto expected = MPMC::PhasePresence(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit);

    require(update.missingPhaseUnstable,
            "SW test stub must certify the missing gas instability");
    require(phaseState.phasePresence.bits() == expected.bits(),
            "SW TPD-triggered reflash must expand O-only to O+G, not unrelated W-only");
    require(primary[Indices::Primary::vaporSaturation] > 0.0,
            "SW TPD-triggered reflash must activate the certified gas phase");
    require(primary[Indices::Primary::waterSaturation] == 0.0,
            "SW TPD-triggered reflash must not manufacture an untriggered water phase");
    require(TriggeredSwFlashStub::restrictedCalls == 1,
            "SW TPD-triggered reflash must first solve exactly one restricted expanded set");
    require(TriggeredSwFlashStub::unrestrictedCalls == 0,
            "stable SW TPD-expanded set must avoid unrelated unrestricted single-role reflash");
}

'''
replace_once(
    test,
    "void checkDependentCompositionCancellationBoundary()\n",
    new_test + "void checkDependentCompositionCancellationBoundary()\n")
replace_once(
    test,
    "        checkTransientNewtonStabilityUsesRestrictedEquilibriumReference();\n        checkDependentCompositionCancellationBoundary();",
    "        checkTransientNewtonStabilityUsesRestrictedEquilibriumReference();\n        checkSwTriggeredReflashPreservesTpdPhaseIdentity();\n        checkDependentCompositionCancellationBoundary();")

print("SW single-phase reflash/canonicalization patch applied")
