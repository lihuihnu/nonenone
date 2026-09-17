from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Read-only probe for the cell-540 O+G -> G -> O+G chatter.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
replace_once(eq, "#include <stdexcept>\n", "#include <stdexcept>\n#include <iostream>\n")
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    struct SwOilRemovalAudit final
    {
        Composition overall{};
        bool reducedConverged{false};
        FlashResult reduced{};
        StabilityResult stability{};
        bool fullConverged{false};
        FlashResult full{};
        bool incipientOilAqueousSupported{false};
    };

    [[nodiscard]] SwOilRemovalAudit auditSwOilRemoval(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        SwOilRemovalAudit audit;
        audit.overall = transitionOverallComposition_(primary, phaseState);
        const auto compositions = normalizedPhaseCompositions_(primary);
        const PhasePresence gasOnly = PhasePresence::gasOnly();
        audit.reduced = flash_.flashRestricted(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            audit.overall,
            gasOnly,
            compositions);
        audit.reducedConverged = audit.reduced.converged;
        if (audit.reducedConverged)
        {
            audit.stability = flash_.stabilityTest(
                primary[Indices::Primary::pressure],
                fluid_.temperature,
                audit.overall,
                audit.reduced.presence,
                audit.reduced.composition);
            audit.incipientOilAqueousSupported =
                fluid_.eos.aqueousVolumeCompositionSupported(
                    audit.stability.incipientComposition[0]);
        }
        audit.full = flash_.flash(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            audit.overall);
        audit.fullConverged = audit.full.converged;
        return audit;
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
                    currentTime_ / 86400.0 >= 0.00550 &&
                    currentTime_ / 86400.0 <= 0.00558 &&
                    pending.phaseState.phasePresence.bits() == 3 &&
                    pending.primary[Indices::Primary::liquidSaturation] <= 0.005;
                if (auditCell540)
                {
                    const auto audit = kernel_.phaseEquilibrium().auditSwOilRemoval(
                        pending.primary, pending.phaseState);
                    std::cerr << "[CELL540-OIL-REMOVAL] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " p_bar=" << pending.primary[Indices::Primary::pressure] / 1.0e5
                              << " So=" << pending.primary[Indices::Primary::liquidSaturation]
                              << " Sg=" << pending.primary[Indices::Primary::vaporSaturation]
                              << " reduced_ok=" << audit.reducedConverged;
                    if (audit.reducedConverged)
                    {
                        std::cerr << " reduced_bits=" << static_cast<int>(audit.reduced.presence.bits())
                                  << " reduced_beta=" << audit.reduced.phaseMoleFraction[0] << ","
                                  << audit.reduced.phaseMoleFraction[1] << ","
                                  << audit.reduced.phaseMoleFraction[2]
                                  << " oil_unstable=" << audit.stability.missingPhaseUnstable[0]
                                  << " oil_sum=" << audit.stability.trialSum[0]
                                  << " oil_excess=" << audit.stability.trialSum[0] - 1.0
                                  << " oil_aq=" << audit.incipientOilAqueousSupported
                                  << " xGred=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << audit.reduced.composition[1][static_cast<std::size_t>(c)];
                        std::cerr << " xOilInc=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << audit.stability.incipientComposition[0][static_cast<std::size_t>(c)];
                    }
                    std::cerr << " full_ok=" << audit.fullConverged;
                    if (audit.fullConverged)
                    {
                        std::cerr << " full_bits=" << static_cast<int>(audit.full.presence.bits())
                                  << " full_S=" << audit.full.saturation[0] << ","
                                  << audit.full.saturation[1] << ","
                                  << audit.full.saturation[2]
                                  << " full_beta=" << audit.full.phaseMoleFraction[0] << ","
                                  << audit.full.phaseMoleFraction[1] << ","
                                  << audit.full.phaseMoleFraction[2];
                    }
                    std::cerr << '\n';
                }

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);

                if (auditCell540)
                {
                    std::cerr << "[CELL540-OIL-REMOVAL-AFTER] after_bits="
                              << static_cast<int>(pending.phaseState.phasePresence.bits())
                              << " status=" << static_cast<int>(update.status)
                              << " removed=" << update.phaseRemoved
                              << " missing_unstable=" << update.missingPhaseUnstable
                              << " So=" << pending.primary[Indices::Primary::liquidSaturation]
                              << " Sg=" << pending.primary[Indices::Primary::vaporSaturation]
                              << " Sw=" << pending.primary[Indices::Primary::waterSaturation]
                              << '\n';
                }
'''
replace_once(state, state_anchor, state_repl)

runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iostream>\n")

print("cell 540 Oil-removal certificate audit injected")
