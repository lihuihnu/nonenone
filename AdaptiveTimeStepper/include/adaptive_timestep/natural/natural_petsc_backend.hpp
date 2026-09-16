/**
 * @file natural_petsc_backend.hpp
 * @brief 连接 AdaptiveTimeStepper 与 Natural/PETSc 运行时的事务后端。
 */
#pragma once

#include <adaptive_timestep/core/solve_result.hpp>
#include <adaptive_timestep/petsc/petsc_error.hpp>
#include <adaptive_timestep/petsc/snes_driver.hpp>
#include <adaptive_timestep/well/control_cycle.hpp>

#include <petscvec.h>

#include <utility>

namespace MPMC
{

struct NoAcceptedStepHook final
{
    template <class Runtime>
    void operator()(Runtime &, Vec, double) const noexcept
    {
    }
};

struct NoFailedSolveHook final
{
    template <class Runtime>
    void operator()(
        Runtime &,
        SNES,
        Vec,
        double,
        double,
        const NonlinearSolveResult &) const noexcept
    {
    }
};

/**
 * @brief Natural PETSc runtime 的事务型自适应时间步后端。
 *
 * 本类只依赖 Runtime 提供的 Natural 公共接口，不依赖 CpGrid 或
 * StructuredGrid。被拒绝的尝试恢复上一接受步的 primary/phase state 与井控；
 * 被接受的尝试通过 runtime.commitTimeStep() 提交历史量。
 *
 * 因此 Land 最大气相饱和度等只在接受步提交的历史不会被失败 Newton 尝试污染。
 */
template <
    class Runtime,
    class WellControlCycleType = NoWellControlCycle,
    class AcceptedStepHook = NoAcceptedStepHook,
    class FailedSolveHook = NoFailedSolveHook>
class NaturalAdaptiveBackend final
{
public:
    NaturalAdaptiveBackend(
        Runtime &runtime,
        SNES snes,
        Vec solution,
        double &currentTime,
        WellControlCycleType wellControlCycle = {},
        AcceptedStepHook acceptedStepHook = {},
        FailedSolveHook failedSolveHook = {},
        bool printNewtonIterations = false,
        NonlinearStagnationConfig stagnationConfig = {})
        : runtime_(runtime),
          solver_(snes, solution, printNewtonIterations, PETSC_COMM_WORLD,
                  stagnationConfig),
          solution_(solution),
          currentTime_(currentTime),
          wellControlCycle_(std::move(wellControlCycle)),
          acceptedStepHook_(std::move(acceptedStepHook)),
          failedSolveHook_(std::move(failedSolveHook))
    {
    }

    [[nodiscard]] double currentTime() const noexcept
    {
        return currentTime_;
    }

    void setTimeStep(double dt)
    {
        runtime_.setTimeStep(dt);
        attemptDt_ = dt;

        // Well schedules are evaluated at the attempted end time.
        runtime_.setCurrentTime(currentTime_ + dt);
    }

    void beginAttempt()
    {
        wellControlCycle_.backup();
    }

    [[nodiscard]] NonlinearSolveResult solve()
    {
        auto result = solver_.solve();
        if (!result.converged)
        {
            // 相出现/消失现在由 Newton 过程中的统一 hysteretic active-set 处理。
            // 这里不再使用“失败后换另一套阈值再重求”的第二套规则；否则一次
            // rescue 删除的相会在下一次 updateState 中按不同门槛重新出现，破坏
            // residual/Jacobian 所对应 active set 的一致性。
            failedSolveHook_(
                runtime_,
                solver_.snes(),
                solution_,
                currentTime_,
                attemptDt_,
                result);
        }
        return result;
    }

    [[nodiscard]] bool updateWellControls()
    {
        return wellControlCycle_.update();
    }

    void rejectAttempt()
    {
        runtime_.rollbackTimeStepAttempt(solution_);
        wellControlCycle_.restore();
    }

    void acceptAttempt(double acceptedTime)
    {
        runtime_.commitTimeStep(solution_);
        acceptedStepHook_(runtime_, solution_, acceptedTime);
        currentTime_ = acceptedTime;
        runtime_.setCurrentTime(acceptedTime);
    }

    void alignCurrentTime(double targetTime)
    {
        currentTime_ = targetTime;
        runtime_.setCurrentTime(targetTime);
    }

private:
    Runtime &runtime_;
    PetscSnesDriver solver_;
    Vec solution_{nullptr};
    double &currentTime_;
    double attemptDt_{0.0};
    WellControlCycleType wellControlCycle_;
    AcceptedStepHook acceptedStepHook_;
    FailedSolveHook failedSolveHook_;
};

} // namespace MPMC
