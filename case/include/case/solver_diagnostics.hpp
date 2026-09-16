/**
 * @file solver_diagnostics.hpp
 * @brief 算例级 Newton/KSP、自适应步长与性能诊断。
 */
#pragma once

#include <adaptive_timestep/core/events.hpp>
#include <adaptive_timestep/core/statistics.hpp>
#include <common/console.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace MPMC::cases
{

/** @brief 求解器与时间步诊断的运行时开关。 */
struct SolverDiagnosticsOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    bool printTimeSteps{true};
    bool writeHistory{true};
    double secondsPerDisplayUnit{86400.0};
    std::string displayTimeUnit{"day"};
};

/**
 * @brief 记录自适应时间步与非线性求解历史的算例级观察器。
 *
 * Live Newton/KSP rows are printed by PetscSnesDriver.  This observer records
 * every nonlinear solve, every Newton row, and every accepted/rejected internal
 * time step in compact CSV files suitable for Python/Matlab post-processing.
 */
class SolverDiagnosticsObserver final
{
public:
    SolverDiagnosticsObserver(PetscMPIInt rank, SolverDiagnosticsOptions options)
        : rank_(rank), options_(std::move(options))
    {
        if (rank_ != 0 || !options_.writeHistory)
            return;

        std::filesystem::create_directories(options_.resultDirectory);
        create_(
            "nonlinear_solve_history.csv",
            "solve_id,start_time,end_time,dt,retry,control_iteration,converged,"
            "snes_iterations,ksp_iterations,initial_residual,final_residual,residual_reduction,"
            "wall_time,reason\n");
        create_(
            "newton_history.csv",
            "solve_id,newton_iteration,residual,residual_ratio,ksp_iterations,ksp_total,"
            "newton_wall_time,ksp_reason\n");
        create_(
            "time_step_history.csv",
            "event_id,status,time,dt,next_dt,retry,snes_iterations,ksp_iterations,solve_wall_time,"
            "well_control_resolves,controls_unsettled,reason\n");
    }

    void fixedTarget(std::size_t step, double targetTime) const noexcept
    {
        if (printEnabled_())
            PetscPrintf(PETSC_COMM_SELF,
                        "\n[TIME][TARGET] output=%zu  target=%.10g %s\n",
                        step,
                        displayTime_(targetTime),
                        options_.displayTimeUnit.c_str());
    }

    void attempt(const MPMC::AdaptiveAttemptEvent &event) const noexcept
    {
        if (printEnabled_() && event.retryCount > 0)
            PetscPrintf(PETSC_COMM_SELF,
                        "[TIME][RETRY ] retry=%d  interval=%.10g -> %.10g %s  dt=%.10g %s\n",
                        event.retryCount,
                        displayTime_(event.startTime),
                        displayTime_(event.endTime),
                        options_.displayTimeUnit.c_str(),
                        displayTime_(event.timeStep),
                        options_.displayTimeUnit.c_str());
    }

    void nonlinearSolve(const MPMC::AdaptiveNonlinearSolveEvent &event) const
    {
        ++solveSerial_;
        if (rank_ != 0 || !options_.writeHistory)
            return;

        const auto &r = event.result;
        {
            auto stream = append_("nonlinear_solve_history.csv");
            stream << std::setprecision(16)
                   << solveSerial_ << ','
                   << displayTime_(event.startTime) << ','
                   << displayTime_(event.endTime) << ','
                   << displayTime_(event.timeStep) << ','
                   << event.retryCount << ','
                   << event.controlIteration << ','
                   << (r.converged ? 1 : 0) << ','
                   << r.nonlinearIterations << ','
                   << r.linearIterations << ','
                   << r.initialResidualNorm << ','
                   << r.finalResidualNorm << ','
                   << r.residualReduction << ','
                   << r.wallTimeSeconds << ','
                   << r.reason << '\n';
        }

        auto newton = append_("newton_history.csv");
        newton << std::setprecision(16);
        for (const auto &row : r.newtonHistory)
        {
            newton << solveSerial_ << ','
                   << row.iteration << ','
                   << row.residualNorm << ','
                   << row.residualRatio << ','
                   << row.linearIterations << ','
                   << row.cumulativeLinearIterations << ','
                   << row.wallTimeSeconds << ','
                   << row.linearReason << '\n';
        }
    }

    void controlSwitch(std::size_t current, std::size_t maximum) const noexcept
    {
        if (printEnabled_())
            PetscPrintf(PETSC_COMM_SELF,
                        "[TIME][WELL  ] control change -> nonlinear re-solve (%zu/%zu)\n",
                        current, maximum);
    }

    void accepted(const MPMC::AdaptiveAcceptedEvent &event) const
    {
        if (printEnabled_())
            PetscPrintf(
                PETSC_COMM_SELF,
                "[TIME][ACCEPT] t=%.10g %s  dt=%.10g %s  SNES=%d  KSP=%lld  solve=%.6f s  well-resolve=%zu  retry=%d\n",
                displayTime_(event.time), options_.displayTimeUnit.c_str(),
                displayTime_(event.timeStep), options_.displayTimeUnit.c_str(),
                event.nonlinearIterations,
                event.linearIterations,
                event.solveTimeSeconds,
                event.controlResolves,
                event.retryCount);

        if (rank_ == 0 && options_.writeHistory)
        {
            ++timeEventSerial_;
            auto stream = append_("time_step_history.csv");
            stream << std::setprecision(16)
                   << timeEventSerial_ << ",ACCEPT,"
                   << displayTime_(event.time) << ','
                   << displayTime_(event.timeStep) << ",,"
                   << event.retryCount << ','
                   << event.nonlinearIterations << ','
                   << event.linearIterations << ','
                   << event.solveTimeSeconds << ','
                   << event.controlResolves << ",0,\n";
        }
    }

    void rejected(const MPMC::AdaptiveRejectedEvent &event) const
    {
        const char *reason = event.result.reason.empty()
            ? "unknown"
            : event.result.reason.c_str();

        if (printEnabled_())
            PetscPrintf(
                PETSC_COMM_SELF,
                "[TIME][REJECT] t=%.10g %s  dt=%.10g -> %.10g %s  retry=%d  SNES=%d  KSP=%lld  solve=%.6f s  reason=%s\n",
                displayTime_(event.time), options_.displayTimeUnit.c_str(),
                displayTime_(event.failedTimeStep),
                displayTime_(event.nextTimeStep), options_.displayTimeUnit.c_str(),
                event.retryCount,
                event.result.nonlinearIterations,
                event.result.linearIterations,
                event.result.wallTimeSeconds,
                reason);

        if (rank_ == 0 && options_.writeHistory)
        {
            ++timeEventSerial_;
            auto stream = append_("time_step_history.csv");
            stream << std::setprecision(16)
                   << timeEventSerial_ << ",REJECT,"
                   << displayTime_(event.time) << ','
                   << displayTime_(event.failedTimeStep) << ','
                   << displayTime_(event.nextTimeStep) << ','
                   << event.retryCount << ','
                   << event.result.nonlinearIterations << ','
                   << event.result.linearIterations << ','
                   << event.result.wallTimeSeconds << ",0,"
                   << (event.controlsUnsettled ? 1 : 0) << ','
                   << reason << '\n';
        }
    }

private:
    [[nodiscard]] bool printEnabled_() const noexcept
    {
        return rank_ == 0 && options_.printTimeSteps;
    }

    [[nodiscard]] double displayTime_(double seconds) const noexcept
    {
        return seconds / options_.secondsPerDisplayUnit;
    }

    void create_(const char *name, const char *header) const
    {
        const auto file = options_.resultDirectory / name;
        std::ofstream stream(file, std::ios::out | std::ios::trunc);
        if (!stream)
            throw std::runtime_error("Failed to create solver diagnostic file: " + file.string());
        stream << header;
    }

    [[nodiscard]] std::ofstream append_(const char *name) const
    {
        const auto file = options_.resultDirectory / name;
        std::ofstream stream(file, std::ios::out | std::ios::app);
        if (!stream)
            throw std::runtime_error("Failed to append solver diagnostic file: " + file.string());
        return stream;
    }

    PetscMPIInt rank_{0};
    SolverDiagnosticsOptions options_;
    mutable std::size_t solveSerial_{0};
    mutable std::size_t timeEventSerial_{0};
};

struct SimulationSummaryInfo final
{
    std::string caseName;
    std::filesystem::path resultDirectory{"./results"};
    std::size_t fixedOutputIntervals{0};
    double finalTimeDays{0.0};
    std::size_t individualWellSwitches{0};
    bool print{true};
    bool write{true};
};

/** @brief 打印并写出便于论文统计的运行结束求解指标。 */
inline void writeSimulationSummary(
    const SimulationSummaryInfo &info,
    const MPMC::AdaptiveTimeStepStatistics &stats,
    double simulationLoopWallSeconds,
    PetscMPIInt rank)
{
    if (rank != 0)
        return;

    PetscMPIInt mpiRanks = 1;
    PetscCallMPIAbort(PETSC_COMM_WORLD, MPI_Comm_size(PETSC_COMM_WORLD, &mpiRanks));

    constexpr double secondsPerDay = 86400.0;
    const double minDtDays = std::isfinite(stats.minimumAcceptedDt)
        ? stats.minimumAcceptedDt / secondsPerDay
        : 0.0;
    const double maxDtDays = stats.maximumAcceptedDt / secondsPerDay;
    const double rejectedSolveTime =
        std::max(0.0, stats.totalSolveTimeSeconds - stats.acceptedSolveTimeSeconds);
    const double rejectedWorkPercent = stats.totalSolveTimeSeconds > 0.0
        ? 100.0 * rejectedSolveTime / stats.totalSolveTimeSeconds
        : 0.0;
    const double avgNPerOutput = info.fixedOutputIntervals == 0
        ? 0.0
        : static_cast<double>(stats.acceptedNonlinearIterations) /
              static_cast<double>(info.fixedOutputIntervals);
    const double avgLPerOutput = info.fixedOutputIntervals == 0
        ? 0.0
        : static_cast<double>(stats.acceptedLinearIterations) /
              static_cast<double>(info.fixedOutputIntervals);
    const double solveShare = simulationLoopWallSeconds > 0.0
        ? 100.0 * stats.totalSolveTimeSeconds / simulationLoopWallSeconds
        : 0.0;

    if (info.print)
    {
        MPMC::ConsoleSection section("SIMULATION SUMMARY");
        section.row("Case", info.caseName)
            .row("MPI ranks", std::to_string(mpiRanks))
            .row("Final time", MPMC::consoleNumber(info.finalTimeDays), "day")
            .row("Fixed output intervals", std::to_string(info.fixedOutputIntervals))
            .separator()
            .row("Accepted / rejected dt",
                 std::to_string(stats.acceptedSteps) + " / " + std::to_string(stats.rejectedSteps))
            .row("Nonlinear solves / re-solves",
                 std::to_string(stats.nonlinearSolves) + " / " + std::to_string(stats.controlSwitches))
            .row("Well-control switches", std::to_string(info.individualWellSwitches))
            .row("Accepted dt [min / max]",
                 MPMC::consoleNumber(minDtDays) + " / " + MPMC::consoleNumber(maxDtDays), "day")
            .separator()
            .row("SNES [attempted / accepted]",
                 std::to_string(stats.attemptedNonlinearIterations) + " / " +
                 std::to_string(stats.acceptedNonlinearIterations))
            .row("KSP [attempted / accepted]",
                 std::to_string(stats.attemptedLinearIterations) + " / " +
                 std::to_string(stats.acceptedLinearIterations))
            .row("Avg SNES / accepted step", MPMC::consoleNumber(stats.averageNonlinearIterationsPerAcceptedStep(), 6))
            .row("Avg KSP / accepted step", MPMC::consoleNumber(stats.averageLinearIterationsPerAcceptedStep(), 6))
            .row("Avg SNES / output", MPMC::consoleNumber(avgNPerOutput, 6))
            .row("Avg KSP / output", MPMC::consoleNumber(avgLPerOutput, 6))
            .row("Avg KSP / Newton", MPMC::consoleNumber(stats.averageLinearIterationsPerNewton(), 6))
            .row("Max SNES / KSP per solve",
                 std::to_string(stats.maximumNonlinearIterationsPerSolve) + " / " +
                 std::to_string(stats.maximumLinearIterationsPerSolve))
            .separator()
            .row("SNES wall [total / accepted]",
                 MPMC::consoleNumber(stats.totalSolveTimeSeconds, 6) + " / " +
                 MPMC::consoleNumber(stats.acceptedSolveTimeSeconds, 6), "s")
            .row("SNES wall [avg / max]",
                 MPMC::consoleNumber(stats.averageSolveTimeSeconds(), 6) + " / " +
                 MPMC::consoleNumber(stats.maximumSolveTimeSeconds, 6), "s")
            .row("Rejected-solve overhead",
                 MPMC::consoleNumber(rejectedSolveTime, 6) + " s (" +
                 MPMC::consoleNumber(rejectedWorkPercent, 3) + "%)")
            .row("Simulation-loop wall", MPMC::consoleNumber(simulationLoopWallSeconds, 6), "s")
            .row("SNES share of loop wall", MPMC::consoleNumber(solveShare, 3), "%");
        const std::string text = section.str();
        PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
    }

    if (!info.write)
        return;

    std::filesystem::create_directories(info.resultDirectory);
    {
        const auto file = info.resultDirectory / "simulation_summary.csv";
        std::ofstream out(file, std::ios::out | std::ios::trunc);
        if (!out)
            throw std::runtime_error("Failed to create simulation_summary.csv.");
        out
            << "case,mpi_ranks,final_time_day,fixed_output_intervals,accepted_internal_steps,rejected_internal_steps,"
            << "nonlinear_solves,well_control_resolves,individual_well_control_switches,min_dt_day,max_dt_day,"
            << "attempted_snes,accepted_snes,attempted_ksp,accepted_ksp,avg_snes_per_internal_step,"
            << "avg_ksp_per_internal_step,avg_snes_per_output,avg_ksp_per_output,avg_ksp_per_newton,"
            << "max_snes_per_solve,max_ksp_per_solve,total_snes_wall,accepted_snes_wall,avg_snes_wall,"
            << "max_snes_wall,rejected_solve_wall,rejected_overhead_percent,simulation_loop_wall,"
            << "snes_solve_share_of_loop_percent\n";
        out << std::setprecision(16)
            << info.caseName << ',' << mpiRanks << ','
            << info.finalTimeDays << ',' << info.fixedOutputIntervals << ','
            << stats.acceptedSteps << ',' << stats.rejectedSteps << ','
            << stats.nonlinearSolves << ',' << stats.controlSwitches << ',' << info.individualWellSwitches << ','
            << minDtDays << ',' << maxDtDays << ','
            << stats.attemptedNonlinearIterations << ',' << stats.acceptedNonlinearIterations << ','
            << stats.attemptedLinearIterations << ',' << stats.acceptedLinearIterations << ','
            << stats.averageNonlinearIterationsPerAcceptedStep() << ','
            << stats.averageLinearIterationsPerAcceptedStep() << ','
            << avgNPerOutput << ',' << avgLPerOutput << ','
            << stats.averageLinearIterationsPerNewton() << ','
            << stats.maximumNonlinearIterationsPerSolve << ',' << stats.maximumLinearIterationsPerSolve << ','
            << stats.totalSolveTimeSeconds << ',' << stats.acceptedSolveTimeSeconds << ','
            << stats.averageSolveTimeSeconds() << ',' << stats.maximumSolveTimeSeconds << ','
            << rejectedSolveTime << ',' << rejectedWorkPercent << ',' << simulationLoopWallSeconds << ','
            << solveShare << '\n';
    }
}

} // namespace MPMC::cases
