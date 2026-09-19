/**
 * @file snes_driver.hpp
 * @brief SNES 求解驱动及 Newton/KSP 性能统计适配器。
 */
#pragma once

#include <adaptive_timestep/core/nonlinear_stagnation.hpp>
#include <adaptive_timestep/core/solve_result.hpp>
#include <adaptive_timestep/petsc/petsc_error.hpp>
#include <common/console.hpp>

#include <petscksp.h>
#include <petscsnes.h>
#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace MPMC
{

/**
 * @brief SNES 轻量驱动器，并统一记录 Newton/KSP 迭代与耗时。
 *
 * The monitor records one row per Newton update:
 *   Newton iteration, nonlinear residual, residual ratio, KSP iterations,
 *   cumulative KSP iterations and wall time of that Newton update.
 *
 * PETSc's built-in -snes_monitor/-ksp_converged_reason remain available, but
 * case/run.sh does not enable them by default because this monitor already
 * presents the information in a single consistent format.
 */
class PetscSnesDriver final
{
public:
    PetscSnesDriver(
        SNES snes,
        Vec solution,
        bool printNewtonIterations = false,
        MPI_Comm communicator = PETSC_COMM_WORLD,
        NonlinearStagnationConfig stagnationConfig = {},
        bool installStagnationConvergenceTest = true)
        : snes_(snes),
          solution_(solution),
          printNewtonIterations_(printNewtonIterations),
          communicator_(communicator),
          stagnationDetector_(stagnationConfig)
    {
        if (snes_ == nullptr || solution_ == nullptr)
            throw std::invalid_argument("PetscSnesDriver requires non-null SNES and solution Vec.");

        throwOnPetscError(
            SNESMonitorSet(snes_, &PetscSnesDriver::monitor_, this, nullptr),
            "SNESMonitorSet(MPMC Newton monitor)");

        if (stagnationDetector_.config().enabled &&
            installStagnationConvergenceTest)
        {
            // Standalone driver compatibility path. NaturalAdaptiveBackend
            // installs one composite convergence callback of its own so the
            // stagnation guard cannot overwrite mesh/mass acceptance gates.
            throwOnPetscError(
                SNESSetConvergenceTest(
                    snes_, &PetscSnesDriver::convergenceTest_, this, nullptr),
                "SNESSetConvergenceTest(MPMC stagnation guard)");
        }
    }

    [[nodiscard]] SNES snes() const noexcept { return snes_; }

    [[nodiscard]] NonlinearSolveResult solve()
    {
        ++solveSerial_;
        history_.clear();
        previousResidualNorm_ = 0.0;
        previousMonitorTime_ = 0.0;
        previousLinearTotal_ = 0;
        monitorInitialized_ = false;
        stagnationDetector_.reset();

        PetscLogDouble start = 0.0;
        PetscLogDouble end = 0.0;
        throwOnPetscError(PetscTime(&start), "PetscTime(start)");

        if (printNewtonIterations_)
        {
            const std::string title =
                consoleCentered("NONLINEAR SOLVE #" + std::to_string(solveSerial_));
            const std::string separator = consoleRule('-');
            throwOnPetscError(
                PetscPrintf(
                    communicator_,
                    "\n%s\n"
                    "  %-7s %-14s %-11s %-7s %-10s %-10s %s\n"
                    "%s\n",
                    title.c_str(),
                    "Newton", "Residual", "Ratio", "KSP", "KSP-total", "Newton(s)", "KSP reason",
                    separator.c_str()),
                "PetscPrintf(SNES begin)");
        }

        throwOnPetscError(SNESSolve(snes_, nullptr, solution_), "SNESSolve");
        throwOnPetscError(PetscTime(&end), "PetscTime(end)");

        SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
        PetscInt nonlinearIterations = 0;
        PetscInt linearIterations = 0;

        throwOnPetscError(SNESGetConvergedReason(snes_, &reason), "SNESGetConvergedReason");
        throwOnPetscError(SNESGetIterationNumber(snes_, &nonlinearIterations), "SNESGetIterationNumber");
        throwOnPetscError(SNESGetLinearSolveIterations(snes_, &linearIterations), "SNESGetLinearSolveIterations");

        NonlinearSolveResult result;
        result.converged = static_cast<int>(reason) > 0;
        result.nonlinearIterations = static_cast<int>(nonlinearIterations);
        result.reasonCode = static_cast<int>(reason);
        result.reason = SNESConvergedReasons[reason];
        const double localWallTime = static_cast<double>(end - start);
        double globalWallTime = localWallTime;
        PetscCallMPIAbort(
            communicator_,
            MPI_Allreduce(&localWallTime, &globalWallTime, 1, MPI_DOUBLE, MPI_MAX, communicator_));
        result.wallTimeSeconds = globalWallTime;
        result.linearIterations = static_cast<long long>(linearIterations);
        result.initialResidualNorm = history_.empty() ? 0.0 : history_.front().residualNorm;
        result.finalResidualNorm = history_.empty() ? 0.0 : history_.back().residualNorm;
        result.residualReduction =
            result.initialResidualNorm > 0.0
                ? result.finalResidualNorm / result.initialResidualNorm
                : 0.0;
        result.newtonHistory = history_;

        if (printNewtonIterations_)
        {
            ConsoleSection section("NONLINEAR SOLVE RESULT #" + std::to_string(solveSerial_));
            section.row("Status", result.converged ? "CONVERGED" : "DIVERGED")
                .row("SNES / KSP iterations",
                     std::to_string(result.nonlinearIterations) + " / " +
                     std::to_string(result.linearIterations))
                .row("Residual r0 -> rf",
                     consoleScientific(result.initialResidualNorm, 6) + " -> " +
                     consoleScientific(result.finalResidualNorm, 6))
                .row("Residual ratio rf/r0", consoleScientific(result.residualReduction, 3))
                .row("Solve wall time", consoleNumber(result.wallTimeSeconds, 6), "s")
                .row("SNES reason", result.reason.empty() ? "unknown" : result.reason);
            const std::string text = section.str();
            throwOnPetscError(
                PetscPrintf(communicator_, "%s", text.c_str()),
                "PetscPrintf(SNES summary)");
        }

        return result;
    }

private:
    static PetscErrorCode convergenceTest_(
        SNES snes,
        PetscInt iteration,
        PetscReal xnorm,
        PetscReal snorm,
        PetscReal fnorm,
        SNESConvergedReason *reason,
        void *context)
    {
        PetscFunctionBegin;

        // PETSc 3.22.2 在 petscsnes.h 中公开 SNESConvergedDefault。必须先调用
        // 默认判据：本 guard 只做补充，不替换用户的 atol/rtol/stol/max-it 设置。
        PetscCall(SNESConvergedDefault(
            snes, iteration, xnorm, snorm, fnorm, reason, nullptr));
        if (*reason != SNES_CONVERGED_ITERATING)
            PetscFunctionReturn(PETSC_SUCCESS);

        auto *self = static_cast<PetscSnesDriver *>(context);
        if (self != nullptr &&
            self->stagnationDetector_.update(
                static_cast<int>(iteration), static_cast<double>(fnorm)))
        {
            *reason = SNES_DIVERGED_LOCAL_MIN;
        }

        PetscFunctionReturn(PETSC_SUCCESS);
    }

    static PetscErrorCode monitor_(
        SNES snes,
        PetscInt iteration,
        PetscReal residualNorm,
        void *context)
    {
        auto *self = static_cast<PetscSnesDriver *>(context);
        return self->recordMonitor_(snes, iteration, residualNorm);
    }

    PetscErrorCode recordMonitor_(
        SNES snes,
        PetscInt iteration,
        PetscReal residualNorm)
    {
        PetscFunctionBegin;

        PetscInt cumulativeLinear = 0;
        PetscLogDouble now = 0.0;
        PetscCall(SNESGetLinearSolveIterations(snes, &cumulativeLinear));
        PetscCall(PetscTime(&now));

        NewtonIterationRecord record;
        record.iteration = static_cast<int>(iteration);
        record.residualNorm = static_cast<double>(residualNorm);
        record.cumulativeLinearIterations = static_cast<long long>(cumulativeLinear);

        if (!monitorInitialized_)
        {
            monitorInitialized_ = true;
            previousResidualNorm_ = record.residualNorm;
            previousMonitorTime_ = static_cast<double>(now);
            previousLinearTotal_ = cumulativeLinear;
            record.residualRatio = 1.0;
            record.linearIterations = 0;
            record.wallTimeSeconds = 0.0;
            record.linearReason = "-";
        }
        else
        {
            record.linearIterations = static_cast<int>(cumulativeLinear - previousLinearTotal_);
            const double localNewtonWall = static_cast<double>(now) - previousMonitorTime_;
            double globalNewtonWall = localNewtonWall;
            PetscCallMPIAbort(
                communicator_,
                MPI_Allreduce(&localNewtonWall, &globalNewtonWall, 1, MPI_DOUBLE, MPI_MAX, communicator_));
            record.wallTimeSeconds = globalNewtonWall;
            record.residualRatio =
                std::abs(previousResidualNorm_) > 0.0
                    ? record.residualNorm / previousResidualNorm_
                    : 0.0;

            KSP ksp = nullptr;
            KSPConvergedReason kspReason = KSP_CONVERGED_ITERATING;
            PetscCall(SNESGetKSP(snes, &ksp));
            if (ksp != nullptr)
            {
                PetscCall(KSPGetConvergedReason(ksp, &kspReason));
                record.linearReasonCode = static_cast<int>(kspReason);
                record.linearReason = KSPConvergedReasons[kspReason];
            }

            previousResidualNorm_ = record.residualNorm;
            previousMonitorTime_ = static_cast<double>(now);
            previousLinearTotal_ = cumulativeLinear;
        }

        history_.push_back(record);

        if (printNewtonIterations_)
        {
            PetscCall(PetscPrintf(
                communicator_,
                "  %-7d % .6e  % .3e  %-7d %-10lld %-10.5f %s\n",
                record.iteration,
                record.residualNorm,
                record.residualRatio,
                record.linearIterations,
                record.cumulativeLinearIterations,
                record.wallTimeSeconds,
                record.linearReason.empty() ? "-" : record.linearReason.c_str()));
        }

        PetscFunctionReturn(PETSC_SUCCESS);
    }

    SNES snes_{nullptr};
    Vec solution_{nullptr};
    bool printNewtonIterations_{false};
    MPI_Comm communicator_{PETSC_COMM_WORLD};

    std::size_t solveSerial_{0};
    std::vector<NewtonIterationRecord> history_;
    double previousResidualNorm_{0.0};
    double previousMonitorTime_{0.0};
    PetscInt previousLinearTotal_{0};
    bool monitorInitialized_{false};
    NonlinearStagnationDetector stagnationDetector_{};
};

} // namespace MPMC
