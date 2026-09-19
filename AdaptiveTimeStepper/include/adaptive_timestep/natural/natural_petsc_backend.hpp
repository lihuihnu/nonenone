/**
 * @file natural_petsc_backend.hpp
 * @brief 连接 AdaptiveTimeStepper 与 Natural/PETSc 运行时的事务后端。
 */
#pragma once

#include <adaptive_timestep/core/solve_result.hpp>
#include <adaptive_timestep/petsc/petsc_error.hpp>
#include <adaptive_timestep/petsc/snes_driver.hpp>
#include <adaptive_timestep/well/control_cycle.hpp>

#include <petscsnes.h>
#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
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
                  stagnationConfig, false),
          solution_(solution),
          currentTime_(currentTime),
          wellControlCycle_(std::move(wellControlCycle)),
          acceptedStepHook_(std::move(acceptedStepHook)),
          failedSolveHook_(std::move(failedSolveHook)),
          stagnationDetector_(stagnationConfig)
    {
        installCompositeConvergenceTest_(snes);
        PetscBool auditAcceptedResidual = PETSC_FALSE;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetBool(
                nullptr, nullptr,
                "-audit_accepted_residual_consistency",
                &auditAcceptedResidual, nullptr));
        auditAcceptedResidualConsistencyEnabled_ =
            auditAcceptedResidual == PETSC_TRUE;

        PetscReal massTolerance = 0.0;
        PetscBool massToleranceSet = PETSC_FALSE;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr,
                "-snes_global_mass_atol",
                &massTolerance, &massToleranceSet));
        if (massToleranceSet == PETSC_TRUE)
            auditedGlobalMassTolerance_ = static_cast<double>(massTolerance);
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
        stagnationDetector_.reset();
        auto result = solver_.solve();

        if (result.converged && auditAcceptedResidualConsistencyEnabled_)
            auditAcceptedResidualConsistency_(result);

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
    static PetscErrorCode compositeConvergenceTest_(
        SNES snes,
        PetscInt iteration,
        PetscReal xNorm,
        PetscReal stepNorm,
        PetscReal functionNorm,
        SNESConvergedReason *reason,
        void *context)
    {
        PetscFunctionBeginUser;
        PetscCheck(
            context != nullptr,
            PetscObjectComm(reinterpret_cast<PetscObject>(snes)),
            PETSC_ERR_ARG_NULL,
            "Natural composite SNES convergence context is null.");

        auto &self =
            *static_cast<NaturalAdaptiveBackend *>(context);

        SNESConvergedReason defaultReason = SNES_CONVERGED_ITERATING;
        PetscCall(SNESConvergedDefault(
            snes,
            iteration,
            xNorm,
            stepNorm,
            functionNorm,
            &defaultReason,
            nullptr));

        // Preserve PETSc hard divergence reasons (NaN, max function count,
        // divergence tolerance, etc.) exactly as reported.
        if (defaultReason < 0)
        {
            *reason = defaultReason;
            PetscFunctionReturn(PETSC_SUCCESS);
        }

        bool acceptedByStandard = defaultReason > 0;

        // Any positive PETSc convergence reason is only provisional when the
        // mesh-independent Natural gates are enabled. Relative convergence is
        // therefore unable to bypass the RMS/Linf/component-mass conditions.
        if (acceptedByStandard && self.meshGateEnabled_)
        {
            const PetscReal rms =
                self.globalEquationCount_ > 0
                    ? functionNorm /
                        std::sqrt(static_cast<PetscReal>(
                            self.globalEquationCount_))
                    : std::numeric_limits<PetscReal>::infinity();

            if (rms > self.rmsAbsoluteTolerance_)
            {
                acceptedByStandard = false;
            }
            else
            {
                Vec residual = nullptr;
                PetscCall(SNESGetFunction(
                    snes, &residual, nullptr, nullptr));
                PetscCheck(
                    residual != nullptr,
                    PetscObjectComm(reinterpret_cast<PetscObject>(snes)),
                    PETSC_ERR_ARG_NULL,
                    "Natural composite convergence requires the SNES residual vector.");

                PetscReal infinityNorm = 0.0;
                PetscCall(VecNorm(
                    residual, NORM_INFINITY, &infinityNorm));

                if (infinityNorm >
                    self.infinityAbsoluteTolerance_)
                {
                    acceptedByStandard = false;
                }
                else
                {
                    const auto massResidual =
                        self.runtime_.evaluateGlobalSignedMassResidual(
                            residual);
                    if (massResidual.maximumAbsolute() >
                        self.globalSignedMassAbsoluteTolerance_)
                    {
                        acceptedByStandard = false;
                    }
                }
            }
        }

        if (acceptedByStandard)
        {
            *reason = defaultReason;
            PetscFunctionReturn(PETSC_SUCCESS);
        }

        // PETSc has not produced an acceptable converged state (either it was
        // still iterating or one of the Natural acceptance gates vetoed it).
        // Only now may the plateau detector request an early retry with a
        // smaller timestep.
        *reason = SNES_CONVERGED_ITERATING;
        if (self.stagnationDetector_.update(
                static_cast<int>(iteration),
                static_cast<double>(functionNorm)))
        {
            *reason = SNES_DIVERGED_LOCAL_MIN;
        }

        PetscFunctionReturn(PETSC_SUCCESS);
    }

    void installCompositeConvergenceTest_(SNES snes)
    {
        PetscReal rmsTolerance = 0.0;
        PetscReal infinityTolerance = 0.0;
        PetscReal massTolerance = 0.0;
        PetscBool rmsSet = PETSC_FALSE;
        PetscBool infinitySet = PETSC_FALSE;
        PetscBool massSet = PETSC_FALSE;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr, "-snes_rms_atol",
                &rmsTolerance, &rmsSet));
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr, "-snes_linf_atol",
                &infinityTolerance, &infinitySet));
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr, "-snes_global_mass_atol",
                &massTolerance, &massSet));

        const int requestedCount =
            (rmsSet == PETSC_TRUE ? 1 : 0) +
            (infinitySet == PETSC_TRUE ? 1 : 0) +
            (massSet == PETSC_TRUE ? 1 : 0);

        if (requestedCount != 0 && requestedCount != 3)
        {
            throw std::invalid_argument(
                "-snes_rms_atol, -snes_linf_atol and "
                "-snes_global_mass_atol must be supplied together.");
        }

        meshGateEnabled_ = requestedCount == 3;
        if (meshGateEnabled_)
        {
            if (!(rmsTolerance > 0.0) ||
                !std::isfinite(rmsTolerance) ||
                !(infinityTolerance > 0.0) ||
                !std::isfinite(infinityTolerance) ||
                !(massTolerance > 0.0) ||
                !std::isfinite(massTolerance))
            {
                throw std::invalid_argument(
                    "Natural RMS/Linf/global-mass convergence tolerances "
                    "must be finite and positive.");
            }

            Vec residual = nullptr;
            PetscCallAbort(
                PETSC_COMM_WORLD,
                SNESGetFunction(
                    snes, &residual, nullptr, nullptr));
            if (residual == nullptr)
                throw std::logic_error(
                    "Natural composite convergence requires a residual vector.");

            PetscCallAbort(
                PETSC_COMM_WORLD,
                VecGetSize(
                    residual, &globalEquationCount_));
            if (globalEquationCount_ <= 0)
                throw std::logic_error(
                    "Natural composite convergence requires nonzero equations.");

            rmsAbsoluteTolerance_ =
                static_cast<double>(rmsTolerance);
            infinityAbsoluteTolerance_ =
                static_cast<double>(infinityTolerance);
            globalSignedMassAbsoluteTolerance_ =
                static_cast<double>(massTolerance);

            PetscReal oldAtol = 0.0;
            PetscReal rtol = 0.0;
            PetscReal stol = 0.0;
            PetscInt maxIterations = 0;
            PetscInt maxFunctions = 0;
            PetscCallAbort(
                PETSC_COMM_WORLD,
                SNESGetTolerances(
                    snes,
                    &oldAtol,
                    &rtol,
                    &stol,
                    &maxIterations,
                    &maxFunctions));

            const PetscReal l2Tolerance =
                rmsTolerance *
                std::sqrt(static_cast<PetscReal>(
                    globalEquationCount_));

            PetscCallAbort(
                PETSC_COMM_WORLD,
                SNESSetTolerances(
                    snes,
                    l2Tolerance,
                    rtol,
                    stol,
                    maxIterations,
                    maxFunctions));
        }

        // Final ownership point: after this call there is exactly one Natural
        // convergence callback. PetscSnesDriver was explicitly constructed
        // with its standalone stagnation callback disabled.
        PetscCallAbort(
            PETSC_COMM_WORLD,
            SNESSetConvergenceTest(
                snes,
                &NaturalAdaptiveBackend::compositeConvergenceTest_,
                this,
                nullptr));

        PetscPrintf(
            PETSC_COMM_WORLD,
            "[SNES][COMPOSITE-CONV] mesh_gates=%s stagnation=%s "
            "RMS_ATOL=%.12e LINF_ATOL=%.12e "
            "GLOBAL_MASS_ATOL=%.12e kg/s N=%lld\n",
            meshGateEnabled_ ? "ON" : "OFF",
            stagnationDetector_.config().enabled ? "ON" : "OFF",
            rmsAbsoluteTolerance_,
            infinityAbsoluteTolerance_,
            globalSignedMassAbsoluteTolerance_,
            static_cast<long long>(globalEquationCount_));
    }

    void auditAcceptedResidualConsistency_(
        const NonlinearSolveResult &result)
    {
        SNES snes = solver_.snes();
        Vec residual = nullptr;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            SNESGetFunction(snes, &residual, nullptr, nullptr));
        if (residual == nullptr)
            throw std::logic_error(
                "Accepted-residual audit requires the SNES residual vector.");

        const auto cachedMass =
            runtime_.evaluateGlobalSignedMassResidual(residual);

        PetscReal cachedL2 = 0.0;
        PetscReal cachedInf = 0.0;
        PetscCallAbort(PETSC_COMM_WORLD, VecNorm(residual, NORM_2, &cachedL2));
        PetscCallAbort(PETSC_COMM_WORLD, VecNorm(residual, NORM_INFINITY, &cachedInf));

        // Re-evaluate F(X) after SNES has returned, using the exact final
        // solution AND current independent phase-state vector. This exposes a
        // stale residual if the final post-check changed phase state after the
        // residual on which convergence was certified.
        PetscCallAbort(
            PETSC_COMM_WORLD,
            runtime_.formFunction(snes, solution_, residual));

        const auto freshMass =
            runtime_.evaluateGlobalSignedMassResidual(residual);

        PetscReal freshL2 = 0.0;
        PetscReal freshInf = 0.0;
        PetscInt globalRows = 0;
        PetscCallAbort(PETSC_COMM_WORLD, VecNorm(residual, NORM_2, &freshL2));
        PetscCallAbort(PETSC_COMM_WORLD, VecNorm(residual, NORM_INFINITY, &freshInf));
        PetscCallAbort(PETSC_COMM_WORLD, VecGetSize(residual, &globalRows));
        const double freshRms =
            globalRows > 0
                ? static_cast<double>(freshL2) /
                    std::sqrt(static_cast<double>(globalRows))
                : std::numeric_limits<double>::infinity();

        const double cachedMassMax = cachedMass.maximumAbsolute();
        const double freshMassMax = freshMass.maximumAbsolute();
        const double massDelta = freshMassMax - cachedMassMax;

        ++acceptedResidualAuditCount_;
        maximumFreshAcceptedMassResidual_ =
            std::max(maximumFreshAcceptedMassResidual_, freshMassMax);
        maximumAcceptedMassResidualIncrease_ =
            std::max(maximumAcceptedMassResidualIncrease_, massDelta);

        const bool violatesMassGate =
            auditedGlobalMassTolerance_ > 0.0 &&
            freshMassMax > auditedGlobalMassTolerance_;
        const bool materiallyChanged =
            freshMassMax >
                std::max(
                    10.0 * std::numeric_limits<double>::epsilon(),
                    cachedMassMax * 1.01 + 1.0e-15) ||
            std::abs(static_cast<double>(freshL2 - cachedL2)) >
                1.0e-12 * std::max(1.0, static_cast<double>(cachedL2));

        if (violatesMassGate || materiallyChanged)
        {
            PetscPrintf(
                PETSC_COMM_WORLD,
                "[ACCEPTED-RESIDUAL-AUDIT] solve=%zu reason=%s "
                "cached_L2=%.12e fresh_L2=%.12e fresh_RMS=%.12e "
                "cached_Linf=%.12e fresh_Linf=%.12e "
                "cached_mass_max=%.12e fresh_mass_max=%.12e "
                "mass_tol=%.12e violation=%d\n",
                acceptedResidualAuditCount_,
                result.reason.c_str(),
                static_cast<double>(cachedL2),
                static_cast<double>(freshL2),
                freshRms,
                static_cast<double>(cachedInf),
                static_cast<double>(freshInf),
                cachedMassMax,
                freshMassMax,
                auditedGlobalMassTolerance_,
                violatesMassGate ? 1 : 0);
        }
    }

    Runtime &runtime_;
    PetscSnesDriver solver_;
    Vec solution_{nullptr};
    double &currentTime_;
    double attemptDt_{0.0};
    WellControlCycleType wellControlCycle_;
    AcceptedStepHook acceptedStepHook_;
    FailedSolveHook failedSolveHook_;
    NonlinearStagnationDetector stagnationDetector_{};
    bool meshGateEnabled_{false};
    double rmsAbsoluteTolerance_{0.0};
    double infinityAbsoluteTolerance_{0.0};
    double globalSignedMassAbsoluteTolerance_{0.0};
    PetscInt globalEquationCount_{0};
    bool auditAcceptedResidualConsistencyEnabled_{false};
    double auditedGlobalMassTolerance_{0.0};
    std::size_t acceptedResidualAuditCount_{0};
    double maximumFreshAcceptedMassResidual_{0.0};
    double maximumAcceptedMassResidualIncrease_{0.0};
};

} // namespace MPMC
