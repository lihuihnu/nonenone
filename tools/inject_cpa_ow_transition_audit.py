from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Add a read-only thermodynamic probe.  It does not alter production decisions.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    struct CpaTransitionAudit final
    {
        StabilityResult stability{};
        FlashResult unrestricted{};
        bool unrestrictedEvaluated{false};
        Composition overall{};
        std::array<Composition, 3> composition{};
    };

    [[nodiscard]] CpaTransitionAudit auditCpaTransition(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        CpaTransitionAudit audit;
        audit.overall = transitionOverallComposition_(primary, phaseState);
        audit.composition = normalizedPhaseCompositions_(primary);
        audit.stability = flash_.stabilityTest(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            audit.overall,
            phaseState.phasePresence,
            audit.composition);

        bool strongMissing = false;
        if (audit.stability.valid)
        {
            for (CompositionalPhase phase : phases_)
            {
                if (phaseState.phasePresence.contains(phase))
                    continue;
                const std::size_t p = static_cast<std::size_t>(phaseIndex(phase));
                const double margin = phaseState.phaseSuppression.contains(phase)
                    ? NaturalNumerics::phaseHysteresisReappearanceMargin
                    : NaturalNumerics::phaseAppearanceStabilityMargin;
                strongMissing = strongMissing ||
                    (audit.stability.missingPhaseUnstable[p] &&
                     audit.stability.trialSum[p] > 1.0 + margin);
            }
        }

        const std::array<double, 3> saturation{
            primary[Indices::Primary::liquidSaturation],
            primary[Indices::Primary::vaporSaturation],
            primary[Indices::Primary::waterSaturation]};
        bool activeNearBoundary = false;
        for (CompositionalPhase phase : phases_)
        {
            if (!phaseState.phasePresence.contains(phase))
                continue;
            activeNearBoundary = activeNearBoundary ||
                saturation[static_cast<std::size_t>(phaseIndex(phase))] <= 1.0e-3;
        }

        if (strongMissing || activeNearBoundary)
        {
            audit.unrestricted = flash_.flash(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                audit.overall);
            audit.unrestrictedEvaluated = true;
        }
        return audit;
    }

'''
replace_once(eq, anchor, probe + anchor)

runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iostream>\n")

state = "models/include/natural/petsc/detail/natural_petsc_runtime_state.inc"
old = r'''                kernel_.phaseEquilibrium().sanitizePrimaryBeforeFlash(
                    pending.primary, options_.useVariableBounds);

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);
'''
new = r'''                kernel_.phaseEquilibrium().sanitizePrimaryBeforeFlash(
                    pending.primary, options_.useVariableBounds);

                const bool auditSelectedCell =
                    static_cast<long long>(cell) == 420 ||
                    static_cast<long long>(cell) == 480 ||
                    static_cast<long long>(cell) == 540;
                const double auditTimeDay = currentTime_ / 86400.0;
                const bool auditWindow =
                    auditTimeDay >= 0.02050 && auditTimeDay <= 0.02077;

                decltype(kernel_.phaseEquilibrium().auditCpaTransition(
                    pending.primary, pending.phaseState)) cpaAudit{};
                const auto beforePrimary = pending.primary;
                const auto beforePhaseState = pending.phaseState;
                if (auditSelectedCell && auditWindow)
                {
                    cpaAudit = kernel_.phaseEquilibrium().auditCpaTransition(
                        pending.primary, pending.phaseState);
                }

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);

                if (auditSelectedCell && auditWindow)
                {
                    const auto &st = cpaAudit.stability;
                    bool rawMissing = false;
                    if (st.valid)
                    {
                        for (int p = 0; p < 3; ++p)
                        {
                            const CompositionalPhase phase = static_cast<CompositionalPhase>(p);
                            if (!beforePhaseState.phasePresence.contains(phase) &&
                                st.missingPhaseUnstable[static_cast<std::size_t>(p)])
                                rawMissing = true;
                        }
                    }
                    const bool phaseChanged =
                        beforePhaseState.phasePresence.bits() != pending.phaseState.phasePresence.bits();
                    const bool suppressionChanged =
                        beforePhaseState.phaseSuppression.bits() != pending.phaseState.phaseSuppression.bits();
                    const bool multiphase = beforePhaseState.phasePresence.count() > 1;

                    if (rawMissing || phaseChanged || suppressionChanged || multiphase ||
                        update.phaseRemoved || update.missingPhaseUnstable ||
                        update.stabilityInvalid || update.recoverableFailure())
                    {
                        std::cerr << "[CPA-OW] t_day=" << auditTimeDay
                                  << " target_day=" << (currentTime_ + options_.timeStep) / 86400.0
                                  << " dt_day=" << options_.timeStep / 86400.0
                                  << " cell=" << static_cast<long long>(cell)
                                  << " before_bits=" << static_cast<int>(beforePhaseState.phasePresence.bits())
                                  << " before_supp=" << static_cast<int>(beforePhaseState.phaseSuppression.bits())
                                  << " before_S="
                                  << beforePrimary[Indices::Primary::liquidSaturation] << ","
                                  << beforePrimary[Indices::Primary::vaporSaturation] << ","
                                  << beforePrimary[Indices::Primary::waterSaturation]
                                  << " after_bits=" << static_cast<int>(pending.phaseState.phasePresence.bits())
                                  << " after_supp=" << static_cast<int>(pending.phaseState.phaseSuppression.bits())
                                  << " after_S="
                                  << pending.primary[Indices::Primary::liquidSaturation] << ","
                                  << pending.primary[Indices::Primary::vaporSaturation] << ","
                                  << pending.primary[Indices::Primary::waterSaturation]
                                  << " status=" << static_cast<int>(update.status)
                                  << " removed=" << update.phaseRemoved
                                  << " miss_unstable=" << update.missingPhaseUnstable
                                  << " stab_invalid=" << update.stabilityInvalid
                                  << " restricted_fail=" << update.restrictedFlashFailed
                                  << " unrestricted_fail=" << update.unrestrictedFlashFailed
                                  << " stab_valid=" << st.valid
                                  << " unstable_OGW="
                                  << st.missingPhaseUnstable[0] << ","
                                  << st.missingPhaseUnstable[1] << ","
                                  << st.missingPhaseUnstable[2]
                                  << " trial_OGW="
                                  << st.trialSum[0] << ","
                                  << st.trialSum[1] << ","
                                  << st.trialSum[2];

                        if (cpaAudit.unrestrictedEvaluated)
                        {
                            std::cerr << " full_conv=" << cpaAudit.unrestricted.converged
                                      << " full_bits=" << static_cast<int>(cpaAudit.unrestricted.presence.bits())
                                      << " full_S="
                                      << cpaAudit.unrestricted.saturation[0] << ","
                                      << cpaAudit.unrestricted.saturation[1] << ","
                                      << cpaAudit.unrestricted.saturation[2]
                                      << " full_beta="
                                      << cpaAudit.unrestricted.phaseMoleFraction[0] << ","
                                      << cpaAudit.unrestricted.phaseMoleFraction[1] << ","
                                      << cpaAudit.unrestricted.phaseMoleFraction[2];
                        }

                        std::cerr << " z=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << cpaAudit.overall[static_cast<std::size_t>(c)];
                        std::cerr << " xO=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << cpaAudit.composition[0][static_cast<std::size_t>(c)];
                        std::cerr << " xW=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << cpaAudit.composition[2][static_cast<std::size_t>(c)];
                        std::cerr << '\n';
                    }
                }
'''
replace_once(state, old, new)
print("CPA O/W transition audit injected")
