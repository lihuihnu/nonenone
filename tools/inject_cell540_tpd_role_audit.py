from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Expose a read-only SW role/stability probe.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
replace_once(eq, "#include <stdexcept>\n", "#include <stdexcept>\n#include <iostream>\n")
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    struct SwTpdRoleAudit final
    {
        StabilityResult stability{};
        Composition overall{};
        double overallOilZ{0.0};
        double overallGasZ{0.0};
        double incipientGasOilZ{0.0};
        double incipientGasGasZ{0.0};
        double gasCompositionDistance{0.0};
        double waterCompositionDistance{0.0};
        bool incipientGasAqueousSupported{false};
        bool incipientWaterAqueousSupported{false};
    };

    [[nodiscard]] SwTpdRoleAudit auditSwTpdRole(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        SwTpdRoleAudit audit;
        audit.overall = transitionOverallComposition_(primary, phaseState);
        const auto compositions = normalizedPhaseCompositions_(primary);
        audit.stability = flash_.stabilityTest(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            audit.overall,
            phaseState.phasePresence,
            compositions);
        if (!fluid_.eos.usesSoreideWhitson())
            return audit;

        const double p = primary[Indices::Primary::pressure];
        const auto overallOil = fluid_.eos.phaseResult(
            p, fluid_.temperature, audit.overall, CompositionalPhase::Oil);
        const auto overallGas = fluid_.eos.phaseResult(
            p, fluid_.temperature, audit.overall, CompositionalPhase::Gas);
        audit.overallOilZ = overallOil.compressibility;
        audit.overallGasZ = overallGas.compressibility;

        const auto &incipientGas = audit.stability.incipientComposition[1];
        const auto incOil = fluid_.eos.phaseResult(
            p, fluid_.temperature, incipientGas, CompositionalPhase::Oil);
        const auto incGas = fluid_.eos.phaseResult(
            p, fluid_.temperature, incipientGas, CompositionalPhase::Gas);
        audit.incipientGasOilZ = incOil.compressibility;
        audit.incipientGasGasZ = incGas.compressibility;
        audit.incipientGasAqueousSupported =
            fluid_.eos.aqueousVolumeCompositionSupported(incipientGas);

        const auto &incipientWater = audit.stability.incipientComposition[2];
        audit.incipientWaterAqueousSupported =
            fluid_.eos.aqueousVolumeCompositionSupported(incipientWater);
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            audit.gasCompositionDistance = std::max(
                audit.gasCompositionDistance,
                std::abs(incipientGas[c] - audit.overall[c]));
            audit.waterCompositionDistance = std::max(
                audit.waterCompositionDistance,
                std::abs(incipientWater[c] - audit.overall[c]));
        }
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
                    currentTime_ / 86400.0 >= 0.00395 &&
                    currentTime_ / 86400.0 <= 0.00405;
                if (auditCell540)
                {
                    const auto audit =
                        kernel_.phaseEquilibrium().auditSwTpdRole(
                            pending.primary, pending.phaseState);
                    const auto &stability = audit.stability;
                    std::cerr << "[CELL540-SWROLE] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " bits=" << static_cast<int>(pending.phaseState.phasePresence.bits())
                              << " p_bar=" << pending.primary[Indices::Primary::pressure] / 1.0e5
                              << " gas_unstable=" << stability.missingPhaseUnstable[1]
                              << " gas_sum=" << stability.trialSum[1]
                              << " gas_excess=" << stability.trialSum[1] - 1.0
                              << " gas_dxinf=" << audit.gasCompositionDistance
                              << " gas_aq_supported=" << audit.incipientGasAqueousSupported
                              << " water_unstable=" << stability.missingPhaseUnstable[2]
                              << " water_sum=" << stability.trialSum[2]
                              << " water_excess=" << stability.trialSum[2] - 1.0
                              << " water_dxinf=" << audit.waterCompositionDistance
                              << " water_aq_supported=" << audit.incipientWaterAqueousSupported
                              << " z_ref=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c == 0 ? "" : ",")
                                  << audit.overall[static_cast<std::size_t>(c)];
                    std::cerr << " gas_x=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c == 0 ? "" : ",")
                                  << stability.incipientComposition[1][static_cast<std::size_t>(c)];
                    std::cerr << " water_x=";
                    for (int c = 0; c < Indices::numComponents; ++c)
                        std::cerr << (c == 0 ? "" : ",")
                                  << stability.incipientComposition[2][static_cast<std::size_t>(c)];
                    std::cerr << " zroot_ref_O/G="
                              << audit.overallOilZ << "/" << audit.overallGasZ
                              << " zroot_inc_O/G="
                              << audit.incipientGasOilZ << "/" << audit.incipientGasGasZ
                              << '\n';
                }

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);

                if (auditCell540)
                {
                    std::cerr << "[CELL540-SWROLE-ACTIVE] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " after_bits=" << static_cast<int>(pending.phaseState.phasePresence.bits())
                              << " status=" << static_cast<int>(update.status)
                              << " missing_unstable=" << update.missingPhaseUnstable
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

print("cell 540 SW TPD role audit injected")
