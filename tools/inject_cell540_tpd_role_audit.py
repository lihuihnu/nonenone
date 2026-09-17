from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Expose a read-only SW role/stability probe around the second cell-540 barrier.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
replace_once(eq, "#include <stdexcept>\n", "#include <stdexcept>\n#include <iostream>\n")
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    struct SwTpdRoleAudit final
    {
        StabilityResult stability{};
        Composition overall{};
        std::array<Composition, 3> compositions{};
        std::array<bool, 3> aqueousSupported{false, false, false};
        double oilGasLogFugacityMismatch{0.0};
        double oilZ{0.0};
        double gasZ{0.0};
        double waterTrialOilZ{0.0};
        double waterTrialGasZ{0.0};
    };

    [[nodiscard]] SwTpdRoleAudit auditSwTpdRole(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        SwTpdRoleAudit audit;
        audit.overall = transitionOverallComposition_(primary, phaseState);
        audit.compositions = normalizedPhaseCompositions_(primary);
        audit.stability = flash_.stabilityTest(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            audit.overall,
            phaseState.phasePresence,
            audit.compositions);
        if (!fluid_.eos.usesSoreideWhitson())
            return audit;

        const double p = primary[Indices::Primary::pressure];
        for (int phase = 0; phase < 3; ++phase)
        {
            audit.aqueousSupported[static_cast<std::size_t>(phase)] =
                fluid_.eos.aqueousVolumeCompositionSupported(
                    audit.compositions[static_cast<std::size_t>(phase)]);
        }

        const auto oil = fluid_.eos.phaseResult(
            p, fluid_.temperature, audit.compositions[0], CompositionalPhase::Oil);
        const auto gas = fluid_.eos.phaseResult(
            p, fluid_.temperature, audit.compositions[1], CompositionalPhase::Gas);
        audit.oilZ = oil.compressibility;
        audit.gasZ = gas.compressibility;
        if (phaseState.phasePresence.contains(CompositionalPhase::Oil) &&
            phaseState.phasePresence.contains(CompositionalPhase::Gas))
        {
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double fo = std::max(oil.fugacity[c], 1.0e-300);
                const double fg = std::max(gas.fugacity[c], 1.0e-300);
                audit.oilGasLogFugacityMismatch = std::max(
                    audit.oilGasLogFugacityMismatch,
                    std::abs(std::log(fo) - std::log(fg)));
            }
        }

        const auto &waterTrial = audit.stability.incipientComposition[2];
        const auto trialOil = fluid_.eos.phaseResult(
            p, fluid_.temperature, waterTrial, CompositionalPhase::Oil);
        const auto trialGas = fluid_.eos.phaseResult(
            p, fluid_.temperature, waterTrial, CompositionalPhase::Gas);
        audit.waterTrialOilZ = trialOil.compressibility;
        audit.waterTrialGasZ = trialGas.compressibility;
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
                    currentTime_ / 86400.0 <= 0.00558;
                if (auditCell540)
                {
                    const auto audit =
                        kernel_.phaseEquilibrium().auditSwTpdRole(
                            pending.primary, pending.phaseState);
                    const auto &stability = audit.stability;
                    std::cerr << "[CELL540-OG-AUDIT] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " bits=" << static_cast<int>(pending.phaseState.phasePresence.bits())
                              << " p_bar=" << pending.primary[Indices::Primary::pressure] / 1.0e5
                              << " So=" << pending.primary[Indices::Primary::liquidSaturation]
                              << " Sg=" << pending.primary[Indices::Primary::vaporSaturation]
                              << " Sw=" << pending.primary[Indices::Primary::waterSaturation]
                              << " og_logf_inf=" << audit.oilGasLogFugacityMismatch
                              << " O_aq=" << audit.aqueousSupported[0]
                              << " G_aq=" << audit.aqueousSupported[1]
                              << " Wslot_aq=" << audit.aqueousSupported[2]
                              << " water_unstable=" << stability.missingPhaseUnstable[2]
                              << " water_sum=" << stability.trialSum[2]
                              << " water_excess=" << stability.trialSum[2] - 1.0
                              << " water_trial_aq="
                              << kernel_.fluid().eos.aqueousVolumeCompositionSupported(
                                     stability.incipientComposition[2])
                              << " zroot_O/G=" << audit.oilZ << "/" << audit.gasZ
                              << " water_trial_zroot_O/G="
                              << audit.waterTrialOilZ << "/" << audit.waterTrialGasZ
                              << " xO=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c == 0 ? "" : ",")
                                  << audit.compositions[0][static_cast<std::size_t>(c)];
                    std::cerr << " xG=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c == 0 ? "" : ",")
                                  << audit.compositions[1][static_cast<std::size_t>(c)];
                    std::cerr << " water_trial_x=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c == 0 ? "" : ",")
                                  << stability.incipientComposition[2][static_cast<std::size_t>(c)];
                    std::cerr << '\n';
                }

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);

                if (auditCell540)
                {
                    std::cerr << "[CELL540-OG-ACTIVE] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " after_bits=" << static_cast<int>(pending.phaseState.phasePresence.bits())
                              << " status=" << static_cast<int>(update.status)
                              << " removed=" << update.phaseRemoved
                              << " missing_unstable=" << update.missingPhaseUnstable
                              << " stability_invalid=" << update.stabilityInvalid
                              << " unrestricted_fail=" << update.unrestrictedFlashFailed
                              << " So=" << pending.primary[Indices::Primary::liquidSaturation]
                              << " Sg=" << pending.primary[Indices::Primary::vaporSaturation]
                              << " Sw=" << pending.primary[Indices::Primary::waterSaturation]
                              << '\n';
                }
'''
replace_once(state, state_anchor, state_repl)

runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iostream>\n")

print("cell 540 O+G SW barrier audit injected")
