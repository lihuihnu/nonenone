/**
 * @file legacy_runtime_initialization.hpp
 * @brief 统一旧油气 formulation 的运行时相态与历史状态初始化流程。
 */
#pragma once

#include <natural/state/phase_state_data.hpp>

namespace MPMC::cases
{

/** @brief 初始化旧 formulation 的 secondary phase state 与 history。 */
template <class Indices, class InitialState, class Runtime, class Solution>
void initializeLegacyRuntime(Runtime &runtime, Solution solution)
{
    const auto initialPhase = MPMC::hydrocarbonPhaseStateFromFlag(
        InitialState::hydrocarbonPhaseFlag);

    if constexpr (InitialState::useLegacySecondaryState)
    {
        MPMC::PhaseStateData<Indices> state;
        state.phase = initialPhase;
        state.equilibriumRatio = InitialState::equilibriumRatio;
        state.overallComposition = InitialState::overallComposition;
        state.liquidMoleFraction = InitialState::liquidMoleFraction;
        runtime.initializeUniformPhaseState(state);
    }
    else
    {
        runtime.initializePhaseStateFromSolution(solution, initialPhase);
        runtime.updateState(solution);
    }
    runtime.initializeHistory(solution);
}

} // namespace MPMC::cases
