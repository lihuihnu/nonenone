/**
 * @file statistics.hpp
 * @brief 自适应时间推进过程的累计求解统计。
 */
#pragma once

#include <adaptive_timestep/core/solve_result.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>

namespace MPMC
{

/**
 * @brief AdaptiveTimeStepper 累积的全局求解工作量统计。
 *
 * "attempted" includes rejected time-step attempts and well-control re-solves;
 * "accepted" counts only nonlinear/linear work belonging to accepted physical
 * time steps.  This distinction is useful when reporting robustness and solver
 * efficiency in papers.
 */
struct AdaptiveTimeStepStatistics final
{
    std::size_t acceptedSteps{0};
    std::size_t rejectedSteps{0};
    std::size_t nonlinearSolves{0};
    // Number of accepted/rejected solve cycles that required a well-control re-solve.
    // Individual well switches are recorded separately by WellControlSwitchOutput.
    std::size_t controlSwitches{0};

    long long attemptedNonlinearIterations{0};
    long long acceptedNonlinearIterations{0};
    long long attemptedLinearIterations{0};
    long long acceptedLinearIterations{0};

    int maximumNonlinearIterationsPerSolve{0};
    long long maximumLinearIterationsPerSolve{0};

    double totalSolveTimeSeconds{0.0};
    double acceptedSolveTimeSeconds{0.0};
    double maximumSolveTimeSeconds{0.0};

    double minimumAcceptedDt{std::numeric_limits<double>::infinity()};
    double maximumAcceptedDt{0.0};

    void recordSolve(const NonlinearSolveResult &result) noexcept
    {
        ++nonlinearSolves;
        attemptedNonlinearIterations += result.nonlinearIterations;
        attemptedLinearIterations += result.linearIterations;
        totalSolveTimeSeconds += result.wallTimeSeconds;
        maximumNonlinearIterationsPerSolve =
            std::max(maximumNonlinearIterationsPerSolve, result.nonlinearIterations);
        maximumLinearIterationsPerSolve =
            std::max(maximumLinearIterationsPerSolve, result.linearIterations);
        maximumSolveTimeSeconds =
            std::max(maximumSolveTimeSeconds, result.wallTimeSeconds);
    }

    void recordAcceptedStep(
        int nonlinearIterations,
        long long linearIterations,
        double solveTimeSeconds,
        double dt) noexcept
    {
        ++acceptedSteps;
        acceptedNonlinearIterations += nonlinearIterations;
        acceptedLinearIterations += linearIterations;
        acceptedSolveTimeSeconds += solveTimeSeconds;
        minimumAcceptedDt = std::min(minimumAcceptedDt, dt);
        maximumAcceptedDt = std::max(maximumAcceptedDt, dt);
    }

    [[nodiscard]] double averageNonlinearIterationsPerAcceptedStep() const noexcept
    {
        return acceptedSteps == 0
            ? 0.0
            : static_cast<double>(acceptedNonlinearIterations) /
                  static_cast<double>(acceptedSteps);
    }

    [[nodiscard]] double averageLinearIterationsPerAcceptedStep() const noexcept
    {
        return acceptedSteps == 0
            ? 0.0
            : static_cast<double>(acceptedLinearIterations) /
                  static_cast<double>(acceptedSteps);
    }

    [[nodiscard]] double averageLinearIterationsPerNewton() const noexcept
    {
        return attemptedNonlinearIterations == 0
            ? 0.0
            : static_cast<double>(attemptedLinearIterations) /
                  static_cast<double>(attemptedNonlinearIterations);
    }

    [[nodiscard]] double averageSolveTimeSeconds() const noexcept
    {
        return nonlinearSolves == 0
            ? 0.0
            : totalSolveTimeSeconds / static_cast<double>(nonlinearSolves);
    }
};

} // namespace MPMC
