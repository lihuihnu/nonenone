/**
 * @file events.hpp
 * @brief 自适应时间步过程的观察事件与回调数据结构。
 */
#pragma once

#include <adaptive_timestep/core/solve_result.hpp>

#include <cstddef>

namespace MPMC
{

struct AdaptiveAttemptEvent final
{
    double startTime{0.0};
    double endTime{0.0};
    double timeStep{0.0};
    int retryCount{0};
    bool clippedToOutputTime{false};
};

/** @brief 每次 SNES 求解后发出，包括井控制切换触发的重复求解。 */
struct AdaptiveNonlinearSolveEvent final
{
    double startTime{0.0};
    double endTime{0.0};
    double timeStep{0.0};
    int retryCount{0};
    std::size_t controlIteration{0};
    NonlinearSolveResult result;
};

struct AdaptiveAcceptedEvent final
{
    double time{0.0};
    double timeStep{0.0};
    int nonlinearIterations{0};
    long long linearIterations{0};
    double solveTimeSeconds{0.0};
    std::size_t controlResolves{0};
    int retryCount{0};
};

struct AdaptiveRejectedEvent final
{
    double time{0.0};
    double failedTimeStep{0.0};
    double nextTimeStep{0.0};
    int retryCount{0};
    NonlinearSolveResult result;
    bool controlsUnsettled{false};
};

struct NullAdaptiveTimeStepObserver final
{
    void fixedTarget(std::size_t, double) const noexcept {}
    void attempt(const AdaptiveAttemptEvent &) const noexcept {}
    void nonlinearSolve(const AdaptiveNonlinearSolveEvent &) const noexcept {}
    void controlSwitch(std::size_t, std::size_t) const noexcept {}
    void accepted(const AdaptiveAcceptedEvent &) const noexcept {}
    void rejected(const AdaptiveRejectedEvent &) const noexcept {}
};

} // namespace MPMC
