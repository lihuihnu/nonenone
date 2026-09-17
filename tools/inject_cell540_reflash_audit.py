from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Expose read-only probes without changing the production thermodynamic path.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
replace_once(eq, "#include <stdexcept>\n", "#include <stdexcept>\n#include <iostream>\n")
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    [[nodiscard]] StabilityResult auditStability(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        const Composition overall = transitionOverallComposition_(primary, phaseState);
        const auto compositions = normalizedPhaseCompositions_(primary);
        return flash_.stabilityTest(
            primary[Indices::Primary::pressure], fluid_.temperature,
            overall, phaseState.phasePresence, compositions);
    }

    [[nodiscard]] Composition auditTransitionOverall(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        return transitionOverallComposition_(primary, phaseState);
    }

    [[nodiscard]] FlashResult auditUnrestrictedFlash(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        return flash_.flash(
            primary[Indices::Primary::pressure], fluid_.temperature,
            transitionOverallComposition_(primary, phaseState));
    }

    [[nodiscard]] FlashResult auditRestrictedFlash(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState,
        PhasePresence allowed) const
    {
        return flash_.flashRestricted(
            primary[Indices::Primary::pressure], fluid_.temperature,
            transitionOverallComposition_(primary, phaseState), allowed);
    }

    [[nodiscard]] StabilityResult auditFlashStability(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState,
        const FlashResult &result) const
    {
        return flash_.stabilityTest(
            primary[Indices::Primary::pressure], fluid_.temperature,
            transitionOverallComposition_(primary, phaseState),
            result.presence, result.composition);
    }

'''
replace_once(eq, anchor, probe + anchor)

state = "models/include/natural/petsc/detail/natural_petsc_runtime_state.inc"
state_anchor = r'''                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);
'''
state_repl = r'''                const bool auditCell540 =
                    static_cast<long long>(cell) == 540 &&
                    currentTime_ / 86400.0 >= 0.00395 &&
                    currentTime_ / 86400.0 <= 0.00402;
                const auto auditBeforeBits = pending.phaseState.phasePresence.bits();
                if (auditCell540)
                {
                    std::cerr << std::setprecision(17);
                    const auto stability = kernel_.phaseEquilibrium().auditStability(
                        pending.primary, pending.phaseState);
                    const double gasMargin =
                        pending.phaseState.phaseSuppression.contains(CompositionalPhase::Gas)
                            ? NaturalNumerics::phaseHysteresisReappearanceMargin
                            : NaturalNumerics::phaseAppearanceStabilityMargin;
                    const bool strongGas =
                        !pending.phaseState.phasePresence.contains(CompositionalPhase::Gas) &&
                        stability.valid && stability.missingPhaseUnstable[1] &&
                        stability.trialSum[1] > 1.0 + gasMargin;

                    std::cerr << "[CELL540-TPD] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " bits=" << static_cast<int>(auditBeforeBits)
                              << " p_bar=" << pending.primary[Indices::Primary::pressure] / 1.0e5
                              << " So=" << pending.primary[Indices::Primary::liquidSaturation]
                              << " Sg=" << pending.primary[Indices::Primary::vaporSaturation]
                              << " Sw=" << pending.primary[Indices::Primary::waterSaturation]
                              << " valid=" << stability.valid
                              << " stable=" << stability.stable
                              << " gas_unstable=" << stability.missingPhaseUnstable[1]
                              << " water_unstable=" << stability.missingPhaseUnstable[2]
                              << " gas_sum=" << stability.trialSum[1]
                              << " water_sum=" << stability.trialSum[2]
                              << " strong_gas=" << strongGas << '\n';

                    if (strongGas)
                    {
                        const auto z = kernel_.phaseEquilibrium().auditTransitionOverall(
                            pending.primary, pending.phaseState);
                        const auto raw = kernel_.phaseEquilibrium().auditUnrestrictedFlash(
                            pending.primary, pending.phaseState);
                        std::cerr << "[CELL540-REFLASH] t_day=" << currentTime_ / 86400.0
                                  << " dt_day=" << options_.timeStep / 86400.0
                                  << " converged=" << raw.converged
                                  << " iterations=" << raw.iterations
                                  << " presence=" << static_cast<int>(raw.presence.bits())
                                  << " z=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",") << z[static_cast<std::size_t>(c)];
                        std::cerr << " beta=" << raw.phaseMoleFraction[0] << ','
                                  << raw.phaseMoleFraction[1] << ',' << raw.phaseMoleFraction[2]
                                  << " sat=" << raw.saturation[0] << ','
                                  << raw.saturation[1] << ',' << raw.saturation[2] << '\n';

                        const std::array<PhasePresence, 3> pairs{{
                            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                                          PhasePresence::bit(CompositionalPhase::Gas)),
                            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                                          PhasePresence::bit(CompositionalPhase::Water)),
                            PhasePresence(PhasePresence::bit(CompositionalPhase::Gas) |
                                          PhasePresence::bit(CompositionalPhase::Water))}};
                        for (const auto allowed : pairs)
                        {
                            const auto restricted =
                                kernel_.phaseEquilibrium().auditRestrictedFlash(
                                    pending.primary, pending.phaseState, allowed);
                            std::cerr << "[CELL540-RESTRICTED] allowed="
                                      << static_cast<int>(allowed.bits())
                                      << " converged=" << restricted.converged
                                      << " presence=" << static_cast<int>(restricted.presence.bits())
                                      << " beta=" << restricted.phaseMoleFraction[0] << ','
                                      << restricted.phaseMoleFraction[1] << ','
                                      << restricted.phaseMoleFraction[2]
                                      << " sat=" << restricted.saturation[0] << ','
                                      << restricted.saturation[1] << ','
                                      << restricted.saturation[2];
                            if (restricted.converged && !restricted.presence.empty())
                            {
                                const auto cert =
                                    kernel_.phaseEquilibrium().auditFlashStability(
                                        pending.primary, pending.phaseState, restricted);
                                std::cerr << " cert_valid=" << cert.valid
                                          << " cert_stable=" << cert.stable
                                          << " trial_sum=" << cert.trialSum[0] << ','
                                          << cert.trialSum[1] << ',' << cert.trialSum[2];
                            }
                            std::cerr << '\n';
                        }
                    }
                }

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);

                if (auditCell540)
                {
                    std::cerr << std::setprecision(17)
                              << "[CELL540-ACTIVE] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " before_bits=" << static_cast<int>(auditBeforeBits)
                              << " after_bits=" << static_cast<int>(pending.phaseState.phasePresence.bits())
                              << " status=" << static_cast<int>(update.status)
                              << " removed=" << update.phaseRemoved
                              << " missing_unstable=" << update.missingPhaseUnstable
                              << " stability_invalid=" << update.stabilityInvalid
                              << " So=" << pending.primary[Indices::Primary::liquidSaturation]
                              << " Sg=" << pending.primary[Indices::Primary::vaporSaturation]
                              << " Sw=" << pending.primary[Indices::Primary::waterSaturation]
                              << '\n';
                }
'''
replace_once(state, state_anchor, state_repl)

runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iostream>\n#include <iomanip>\n")

print("cell 540 high-precision raw/restricted flash audit instrumentation injected")
