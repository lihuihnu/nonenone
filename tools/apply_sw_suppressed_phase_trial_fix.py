#!/usr/bin/env python3
from pathlib import Path


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one anchor, found {count}")
    path.write_text(text.replace(old, new, 1))

root = Path(__file__).resolve().parents[1]
eq = root / "models/include/natural/state/three_phase_equilibrium.hpp"
test = root / "test/src/unit/sw_flash_recovery_test.cpp"

old_block = r'''                    else
                    {
                        unstableMissingPhase =
                            hasStrongMissingPhaseInstability_(
                                stability, reduced.presence,
                                phaseState.phaseSuppression);
                        updateResult.missingPhaseUnstable = unstableMissingPhase;
                    }
                }

                if (!unstableMissingPhase)
'''
new_block = r'''                    else
                    {
                        unstableMissingPhase =
                            hasStrongMissingPhaseInstability_(
                                stability, reduced.presence,
                                phaseState.phaseSuppression);
                        updateResult.missingPhaseUnstable = unstableMissingPhase;

                        // If the phase removed by this very trial is already
                        // strongly unstable in the certified reduced
                        // equilibrium, the proposed active-set deletion is
                        // thermodynamically inadmissible.  Do not project the
                        // same line-search candidate through an unrestricted
                        // flash, which can reinsert the phase at finite
                        // saturation and create O+G -> G -> O+G chatter.
                        // Restore the transaction and let SNES shorten the
                        // Newton step while approaching the boundary from the
                        // active side.  Instability of a different, previously
                        // missing phase still follows the unrestricted-flash
                        // appearance path below.
                        if (hasStrongSuppressedPhaseInstability_(
                                stability, reduced.presence,
                                phaseState.phaseSuppression))
                        {
                            primary = originalPrimary;
                            phaseState = originalPhaseState;
                            updateResult.phaseRemoved = true;
                            updateResult.status =
                                PhaseUpdateStatus::RecoverableThermodynamicFailure;
                            return updateResult;
                        }
                    }
                }

                if (!unstableMissingPhase)
'''
replace_once(eq, old_block, new_block)

helper_anchor = r'''    /** @brief 缺失相只有明显越过稳定性边界时才重新生成，避免边界 active-set 抖动。 */
    [[nodiscard]] static bool hasStrongMissingPhaseInstability_(
'''
helper = r'''    /**
     * @brief 判定本次 active-set 删除的相是否仍被强 TPD 证明必须存在。
     *
     * A suppressed phase is different from an ordinary missing phase: it was
     * active at the start of this update and was removed only because the
     * Newton trial crossed the saturation probe boundary.  If the certified
     * reduced equilibrium immediately says that same phase is strongly
     * unstable, the deletion trial itself is invalid and must be rejected
     * transactionally rather than followed by a finite-saturation reflash.
     */
    [[nodiscard]] static bool hasStrongSuppressedPhaseInstability_(
        const StabilityResult &stability,
        PhasePresence present,
        PhasePresence suppression)
    {
        if (!stability.valid)
            return false;

        for (CompositionalPhase phase : phases_)
        {
            if (present.contains(phase) || !suppression.contains(phase))
                continue;
            const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
            if (stability.missingPhaseUnstable[p] &&
                stability.trialSum[p] >
                    1.0 + NaturalNumerics::phaseHysteresisReappearanceMargin)
            {
                return true;
            }
        }
        return false;
    }

'''
replace_once(eq, helper_anchor, helper + helper_anchor)

replace_once(
    test,
    "#include <natural/compositional_mixture.hpp>\n",
    "#include <natural/compositional_mixture.hpp>\n"
    "#include <natural/fluid_system.hpp>\n"
    "#include <natural/state/three_phase_equilibrium.hpp>\n")
replace_once(
    test,
    "using Flash = MPMC::CubicThreePhaseFlash<Indices>;\n",
    "using Flash = MPMC::CubicThreePhaseFlash<Indices>;\n"
    "using Equilibrium = MPMC::FullyCompositionalThreePhaseEquilibrium<Indices>;\n")

insert_anchor = r'''

} // namespace

int main()
'''
regression = r'''
void checkLmhRejectedOilDeletionRemainsTransactional()
{
    // Real cell-540 O+G line-search state at the second LMH SW barrier.  Oil
    // has crossed the active-set probe, but the certified Gas-only state has
    // a strong missing-Oil TPD signal (trialSum about 1.0088).  Before this
    // fix updatePhaseState() performed a full flash in the same post-check and
    // reinserted Oil at about 3 mol/vol %, creating a discontinuous
    // O+G -> G -> O+G projection.
    auto eos = makeLmhBoundarySw();
    MPMC::FluidSystem<Indices> fluid(
        {700.0, 220.0, 363.0},
        {2.0e-4, 2.0e-5, 5.0e-5},
        std::move(eos),
        {"H2O", "Light", "Middle", "Heavy"},
        {"Oil", "Gas", "Water"},
        653.2);
    fluid.configureFullyCompositionalThreePhase(0);
    const Equilibrium equilibrium(fluid);

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 279.157e5;
    primary[Indices::Primary::liquidSaturation] = -0.00178598;
    primary[Indices::Primary::vaporSaturation] = 1.00178598;
    primary[Indices::Primary::waterSaturation] = 0.0;

    constexpr Composition oil{
        0.743781, 0.0404795, 0.0996408, 0.1160987};
    constexpr Composition gas{
        0.963397, 0.0170619, 0.0136560, 0.0058851};
    constexpr Composition water{1.0, 0.0, 0.0, 0.0};
    for (int c = 0; c < Indices::numIndependentCompositionsPerPhase; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(c);
        primary[static_cast<std::size_t>(Indices::Primary::liquidComposition[i])] = oil[i];
        primary[static_cast<std::size_t>(Indices::Primary::vaporComposition[i])] = gas[i];
        primary[static_cast<std::size_t>(Indices::Primary::waterComposition[i])] = water[i];
    }

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit);
    phaseState.phaseMoleFraction = {0.0168417, 0.9831583, 0.0};
    for (std::size_t c = 0; c < 4; ++c)
    {
        phaseState.overallComposition[c] =
            phaseState.phaseMoleFraction[0] * oil[c] +
            phaseState.phaseMoleFraction[1] * gas[c];
    }

    const auto originalPrimary = primary;
    const auto result = equilibrium.updatePhaseState(primary, phaseState);

    require(result.phaseRemoved,
            "LMH Oil-boundary trial must be recognized as a deletion probe");
    require(result.missingPhaseUnstable,
            "LMH Gas-only certificate must report strongly unstable missing Oil");
    require(result.status == MPMC::PhaseUpdateStatus::RecoverableThermodynamicFailure,
            "strongly unstable just-deleted Oil must reject the line-search trial");
    require(phaseState.phasePresence.bits() ==
                (MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit),
            "rejected Oil deletion must preserve original O+G presence");
    for (std::size_t i = 0; i < primary.size(); ++i)
        require(primary[i] == originalPrimary[i],
                "rejected Oil deletion must restore the original primary transaction");
}
'''
replace_once(test, insert_anchor, "\n" + regression + insert_anchor)
replace_once(
    test,
    "        checkLmhSinglePhaseRejectsWrongRoleGasTpdBasin();\n",
    "        checkLmhSinglePhaseRejectsWrongRoleGasTpdBasin();\n"
    "        checkLmhRejectedOilDeletionRemainsTransactional();\n")

print("suppressed-phase deletion trial rejection patch applied")
