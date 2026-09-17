from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Read-only diagnostic hook inside the equilibrium class so private flash helpers
# are available without changing any production decision.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    struct CpaGasCandidateAudit final
    {
        StabilityResult stability{};
        std::array<Composition, 3> composition{};
        double pressure{0.0};
        double gasDistanceOil{0.0};
        double gasDistanceWater{0.0};
        double gasZOil{0.0};
        double gasZGas{0.0};
        double gasZWater{0.0};
        double gasGOil{0.0};
        double gasGGas{0.0};
        double gasGWater{0.0};
    };

    [[nodiscard]] CpaGasCandidateAudit auditCpaGasCandidate(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        CpaGasCandidateAudit audit;
        audit.pressure = primary[Indices::Primary::pressure];
        audit.composition = normalizedPhaseCompositions_(primary);
        const auto overall = transitionOverallComposition_(primary, phaseState);
        audit.stability = flash_.stabilityTest(
            audit.pressure, fluid_.temperature, overall,
            phaseState.phasePresence, audit.composition);

        const auto &trial = audit.stability.incipientComposition[1];
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            audit.gasDistanceOil = std::max(
                audit.gasDistanceOil, std::abs(trial[c] - audit.composition[0][c]));
            audit.gasDistanceWater = std::max(
                audit.gasDistanceWater, std::abs(trial[c] - audit.composition[2][c]));
        }

        const auto oil = fluid_.eos.phaseResult(
            audit.pressure, fluid_.temperature, trial, CompositionalPhase::Oil);
        const auto gas = fluid_.eos.phaseResult(
            audit.pressure, fluid_.temperature, trial, CompositionalPhase::Gas);
        const auto water = fluid_.eos.phaseResult(
            audit.pressure, fluid_.temperature, trial, CompositionalPhase::Water);
        audit.gasZOil = oil.compressibility;
        audit.gasZGas = gas.compressibility;
        audit.gasZWater = water.compressibility;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double x = std::max(trial[c], 1.0e-300);
            audit.gasGOil += trial[c] *
                (std::log(x) + std::log(std::max(oil.fugacityCoefficient[c], 1.0e-300)));
            audit.gasGGas += trial[c] *
                (std::log(x) + std::log(std::max(gas.fugacityCoefficient[c], 1.0e-300)));
            audit.gasGWater += trial[c] *
                (std::log(x) + std::log(std::max(water.fugacityCoefficient[c], 1.0e-300)));
        }
        return audit;
    }

'''
replace_once(eq, anchor, probe + anchor)

runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iostream>\n")

state = "models/include/natural/petsc/detail/natural_petsc_runtime_state.inc"
old = r'''                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);
'''
new = r'''                if (static_cast<long long>(cell) == 540 &&
                    currentTime_ / 86400.0 >= 0.02075 &&
                    currentTime_ / 86400.0 <= 0.02077 &&
                    pending.phaseState.phasePresence.bits() == 5)
                {
                    const auto audit = kernel_.phaseEquilibrium().auditCpaGasCandidate(
                        pending.primary, pending.phaseState);
                    const auto &st = audit.stability;
                    if (st.valid && st.missingPhaseUnstable[1])
                    {
                        std::cerr << "[CPA-GAS-CAND] t_day=" << currentTime_ / 86400.0
                                  << " dt_day=" << options_.timeStep / 86400.0
                                  << " p_bar=" << audit.pressure / 1.0e5
                                  << " S="
                                  << pending.primary[Indices::Primary::liquidSaturation] << ","
                                  << pending.primary[Indices::Primary::vaporSaturation] << ","
                                  << pending.primary[Indices::Primary::waterSaturation]
                                  << " trial_sum=" << st.trialSum[1]
                                  << " dist_O=" << audit.gasDistanceOil
                                  << " dist_W=" << audit.gasDistanceWater
                                  << " Z_OGW=" << audit.gasZOil << ","
                                  << audit.gasZGas << "," << audit.gasZWater
                                  << " G_OGW=" << audit.gasGOil << ","
                                  << audit.gasGGas << "," << audit.gasGWater
                                  << " xTrial=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << st.incipientComposition[1][static_cast<std::size_t>(c)];
                        std::cerr << " xO=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << audit.composition[0][static_cast<std::size_t>(c)];
                        std::cerr << " xW=";
                        for (int c = 0; c < Indices::numComponents; ++c)
                            std::cerr << (c == 0 ? "" : ",")
                                      << audit.composition[2][static_cast<std::size_t>(c)];
                        std::cerr << '\n';
                    }
                }

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);
'''
replace_once(state, old, new)
print("CPA Gas candidate identity audit injected")
