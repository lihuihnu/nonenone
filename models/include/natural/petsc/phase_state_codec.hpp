/**
 * @file phase_state_codec.hpp
 * @brief 相状态向量与离散 phase-presence 编码之间的转换。
 */
#pragma once

#include <natural/state/phase_state_data.hpp>

#include <array>
#include <cstddef>

namespace MPMC
{

/**
 * @brief PhaseStateData 与 PETSc 每单元连续数组之间的唯一转换入口。
 *
 * Legacy and fully-compositional modes intentionally use different Vec layouts;
 * this codec is the only place that knows those layouts.
 */
template <class Indices>
struct PetscPhaseStateCodec final
{
    using Data = PhaseStateData<Indices>;
    using Array = std::array<double, Indices::numPhaseStateVariables>;

    [[nodiscard]] static Data decode(const Array &values)
    {
        Data state;
        if constexpr (!Indices::fullyCompositionalThreePhase)
        {
            state.phase = decodeHydrocarbonPhaseState(
                values[static_cast<std::size_t>(Indices::PhaseState::flag)]);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const auto c = static_cast<std::size_t>(component);
                state.equilibriumRatio[c] = values[static_cast<std::size_t>(
                    Indices::PhaseState::equilibriumRatio[c])];
                state.overallComposition[c] = values[static_cast<std::size_t>(
                    Indices::PhaseState::overallComposition[c])];
            }
            state.liquidMoleFraction = values[static_cast<std::size_t>(
                Indices::PhaseState::liquidFraction)];
            state.liquidCompressibility = values[static_cast<std::size_t>(
                Indices::PhaseState::liquidCompressibilityFactor)];
            state.vaporCompressibility = values[static_cast<std::size_t>(
                Indices::PhaseState::vaporCompressibilityFactor)];
        }
        else
        {
            state.phasePresence = decodePhasePresence(
                values[static_cast<std::size_t>(Indices::PhaseState::flag)]);
            state.phaseSuppression = PhasePresence(static_cast<std::uint8_t>(
                values[static_cast<std::size_t>(Indices::PhaseState::phaseSuppressionFlag)] + 0.5));
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const auto c = static_cast<std::size_t>(component);
                state.vaporOilEquilibriumRatio[c] = values[static_cast<std::size_t>(
                    Indices::PhaseState::vaporOilEquilibriumRatio[c])];
                state.waterOilEquilibriumRatio[c] = values[static_cast<std::size_t>(
                    Indices::PhaseState::waterOilEquilibriumRatio[c])];
                state.overallComposition[c] = values[static_cast<std::size_t>(
                    Indices::PhaseState::overallComposition[c])];
            }
            for (int p = 0; p < 3; ++p)
            {
                state.phaseMoleFraction[static_cast<std::size_t>(p)] = values[static_cast<std::size_t>(
                    Indices::PhaseState::phaseMoleFraction[static_cast<std::size_t>(p)])];
            }
            state.phaseCompressibility[0] = values[static_cast<std::size_t>(
                Indices::PhaseState::liquidCompressibilityFactor)];
            state.phaseCompressibility[1] = values[static_cast<std::size_t>(
                Indices::PhaseState::vaporCompressibilityFactor)];
            state.phaseCompressibility[2] = values[static_cast<std::size_t>(
                Indices::PhaseState::waterCompressibilityFactor)];
        }
        return state;
    }

    [[nodiscard]] static Array encode(const Data &state)
    {
        Array values{};
        if constexpr (!Indices::fullyCompositionalThreePhase)
        {
            values[static_cast<std::size_t>(Indices::PhaseState::flag)] =
                encodeHydrocarbonPhaseState(state.phase);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const auto c = static_cast<std::size_t>(component);
                values[static_cast<std::size_t>(Indices::PhaseState::equilibriumRatio[c])] =
                    state.equilibriumRatio[c];
                values[static_cast<std::size_t>(Indices::PhaseState::overallComposition[c])] =
                    state.overallComposition[c];
            }
            values[static_cast<std::size_t>(Indices::PhaseState::liquidFraction)] =
                state.liquidMoleFraction;
            values[static_cast<std::size_t>(Indices::PhaseState::liquidCompressibilityFactor)] =
                state.liquidCompressibility;
            values[static_cast<std::size_t>(Indices::PhaseState::vaporCompressibilityFactor)] =
                state.vaporCompressibility;
        }
        else
        {
            values[static_cast<std::size_t>(Indices::PhaseState::flag)] =
                encodePhasePresence(state.phasePresence);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const auto c = static_cast<std::size_t>(component);
                values[static_cast<std::size_t>(Indices::PhaseState::vaporOilEquilibriumRatio[c])] =
                    state.vaporOilEquilibriumRatio[c];
                values[static_cast<std::size_t>(Indices::PhaseState::waterOilEquilibriumRatio[c])] =
                    state.waterOilEquilibriumRatio[c];
                values[static_cast<std::size_t>(Indices::PhaseState::overallComposition[c])] =
                    state.overallComposition[c];
            }
            for (int p = 0; p < 3; ++p)
            {
                values[static_cast<std::size_t>(
                    Indices::PhaseState::phaseMoleFraction[static_cast<std::size_t>(p)])] =
                    state.phaseMoleFraction[static_cast<std::size_t>(p)];
            }
            values[static_cast<std::size_t>(Indices::PhaseState::liquidCompressibilityFactor)] =
                state.phaseCompressibility[0];
            values[static_cast<std::size_t>(Indices::PhaseState::vaporCompressibilityFactor)] =
                state.phaseCompressibility[1];
            values[static_cast<std::size_t>(Indices::PhaseState::waterCompressibilityFactor)] =
                state.phaseCompressibility[2];
            values[static_cast<std::size_t>(Indices::PhaseState::phaseSuppressionFlag)] =
                static_cast<double>(state.phaseSuppression.bits());
        }
        return values;
    }
};

} // namespace MPMC
