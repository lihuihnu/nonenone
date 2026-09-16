/**
 * @file phase_state.hpp
 * @brief 兼容旧油气模型的运行时相状态表示。
 */
#pragma once

#include <natural/thermo/phase_behavior.hpp>

#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 旧油/气 formulation 的运行时相态。
 *
 * Kept unchanged for source/binary-layout compatibility of all existing cases.
 */
enum class HydrocarbonPhaseState : int
{
    TwoPhase = 0,
    LiquidOnly = 1,
    VaporOnly = 2
};

[[nodiscard]] constexpr HydrocarbonPhaseState hydrocarbonPhaseStateFromFlag(int flag)
{
    switch (flag)
    {
    case 1: return HydrocarbonPhaseState::LiquidOnly;
    case 2: return HydrocarbonPhaseState::VaporOnly;
    default: return HydrocarbonPhaseState::TwoPhase;
    }
}

[[nodiscard]] constexpr int toPhaseStateFlag(HydrocarbonPhaseState state) noexcept
{
    return static_cast<int>(state);
}

[[nodiscard]] inline HydrocarbonPhaseState decodeHydrocarbonPhaseState(double flag)
{
    return hydrocarbonPhaseStateFromFlag(static_cast<int>(std::lround(flag)));
}

[[nodiscard]] inline double encodeHydrocarbonPhaseState(HydrocarbonPhaseState state) noexcept
{
    return static_cast<double>(toPhaseStateFlag(state));
}

} // namespace MPMC
