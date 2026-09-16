/**
 * @file phase_update_result.hpp
 * @brief 相态更新的结构化结果，用于把可恢复热力学失败与普通状态更新分离。
 */
#pragma once

namespace MPMC
{

enum class PhaseUpdateStatus
{
    Updated,
    StableReducedSet,
    RecoverableThermodynamicFailure
};

struct PhaseUpdateResult final
{
    PhaseUpdateStatus status{PhaseUpdateStatus::Updated};
    bool phaseRemoved{false};
    bool missingPhaseUnstable{false};
    bool stabilityInvalid{false};
    bool restrictedFlashFailed{false};
    bool unrestrictedFlashFailed{false};

    [[nodiscard]] bool succeeded() const noexcept
    {
        return status != PhaseUpdateStatus::RecoverableThermodynamicFailure;
    }

    [[nodiscard]] bool recoverableFailure() const noexcept
    {
        return status == PhaseUpdateStatus::RecoverableThermodynamicFailure;
    }
};

} // namespace MPMC
