from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one anchor, found {count}")
    p.write_text(text.replace(old, new, 1))

# Expose a read-only stability probe through the existing equilibrium object.
eq = "models/include/natural/state/three_phase_equilibrium.hpp"
replace_once(eq, "#include <stdexcept>\n", "#include <stdexcept>\n#include <iostream>\n")
anchor = "    /** @brief 将独立 P–T–z flash 结果写入 Natural 主变量和 secondary phase-state。 */\n"
probe = r'''    [[nodiscard]] StabilityResult auditStability(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        const Composition overall =
            transitionOverallComposition_(primary, phaseState);
        const auto compositions = normalizedPhaseCompositions_(primary);
        return flash_.stabilityTest(
            primary[Indices::Primary::pressure],
            fluid_.temperature,
            overall,
            phaseState.phasePresence,
            compositions);
    }

'''
replace_once(eq, anchor, probe + anchor)

# Trace only cell 540 and only the narrow physical-time window around 0.004 day.
state = "models/include/natural/petsc/detail/natural_petsc_runtime_state.inc"
state_anchor = r'''                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);
'''
state_repl = r'''                const bool auditCell540 =
                    static_cast<long long>(cell) == 540 &&
                    currentTime_ / 86400.0 >= 0.00380 &&
                    currentTime_ / 86400.0 <= 0.00410;
                const auto auditBeforeBits = pending.phaseState.phasePresence.bits();
                if (auditCell540)
                {
                    const auto stability =
                        kernel_.phaseEquilibrium().auditStability(
                            pending.primary, pending.phaseState);
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
                              << " water_x=";
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
                    std::cerr << "[CELL540-ACTIVE] t_day=" << currentTime_ / 86400.0
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

# Trace the actual RATE-control AD row produced by the injector perforation.
jac = "models/include/natural/petsc/detail/natural_petsc_runtime_jacobian.inc"
jac_anchor = r'''                const auto &result =
                    evaluatePerforation_(
                        well,
                        perforation,
                        bhp[wellIndex],
                        cache[static_cast<std::size_t>(localBlock)]);
'''
jac_repl = jac_anchor + r'''
                if (static_cast<long long>(perforation.currentCellId) == 540 &&
                    well.isRateControlled() &&
                    currentTime_ / 86400.0 >= 0.00380 &&
                    currentTime_ / 86400.0 <= 0.00410)
                {
                    const auto rateEquation = selectControlledRate<Indices>(
                        well.control,
                        result.surfacePhaseRate,
                        result.reservoirPhaseRate);
                    const double rowScale = wellControlEquationScaleFactor_(wellIndex);
                    const auto &auditState = cache[static_cast<std::size_t>(localBlock)].state;
                    std::cerr << "[CELL540-WELLJAC] t_day=" << currentTime_ / 86400.0
                              << " dt_day=" << options_.timeStep / 86400.0
                              << " bits=" << static_cast<int>(auditState.phasePresence.bits())
                              << " rate=" << scalarValue(rateEquation)
                              << " target=" << well.signedTarget()
                              << " p_bar=" << scalarValue(auditState.pressure) / 1.0e5
                              << " bhp_bar=" << scalarValue(bhp[wellIndex]) / 1.0e5
                              << " d_scaled=";
                    for (PetscInt col = 0; col < N; ++col)
                    {
                        const double scaled = rateEquation.derivative(col) *
                            rowScale * variableScaleFactor_(col);
                        std::cerr << (col == 0 ? "" : ",") << scaled;
                    }
                    std::cerr << '\n';
                }
'''
replace_once(jac, jac_anchor, jac_repl)

# NaturalPetscRuntime owns the injected iostream use in its .inc fragments.
runtime = "models/include/natural/petsc/natural_petsc_runtime.hpp"
replace_once(runtime, "#include <functional>\n", "#include <functional>\n#include <iostream>\n")

print("cell 540 audit instrumentation injected")
