/**
 * @file stepper.hpp
 * @brief 与网格和求解器解耦的事务型自适应时间推进器。
 */
#pragma once

#include <adaptive_timestep/core/config.hpp>
#include <adaptive_timestep/core/events.hpp>
#include <adaptive_timestep/core/exceptions.hpp>
#include <adaptive_timestep/core/policy.hpp>
#include <adaptive_timestep/core/statistics.hpp>

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>

namespace MPMC
{

/**
 * @brief 与网格和非线性求解器解耦的事务型自适应时间推进。
 *
 * 后端只需满足以下最小契约：读取当前时间、设置 dt、开始一次尝试、执行非线性求解、
 * 更新井控制，以及接受/拒绝本次尝试。
 *
 * `beginAttempt()` 与 `rejectAttempt()` 构成事务边界：失败尝试中修改的解、相态、井状态
 * 必须能够完全回滚；只有 `acceptAttempt()` 可以提交新的历史量。
 */
template <class Backend, class Observer = NullAdaptiveTimeStepObserver>
class AdaptiveTimeStepper final
{
public:
    AdaptiveTimeStepper(
        Backend &backend,
        AdaptiveTimeStepConfig config,
        Observer observer = {})
        : backend_(backend),
          policy_(std::move(config)),
          observer_(std::move(observer))
    {
    }

    /** @brief 返回已尝试和已接受内部时间步的累计工作量。 */
    [[nodiscard]] const AdaptiveTimeStepStatistics &statistics() const noexcept
    {
        return statistics_;
    }

    /**
     * @brief 推进 `numberOfFixedSteps` 个固定输出区间。
     *
     * 内部自适应步长可以小于 `fixedOutputDt`，但每个输出时刻都精确对齐；
     * `outputCallback` 只观察已经接受并提交的状态。
     */
    template <class OutputCallback>
    void run(std::size_t numberOfFixedSteps, OutputCallback &&outputCallback)
    {
        (void)runUntil(
            numberOfFixedSteps,
            [](std::size_t, double) { return false; },
            std::forward<OutputCallback>(outputCallback));
    }

    /**
     * @brief Advance fixed output intervals until a committed-state predicate is met.
     *
     * The stop predicate is evaluated only after advanceTo_() has accepted all
     * internal steps to the fixed target and after outputCallback has observed
     * the committed state.  This is suitable for cumulative-volume/PVI stops:
     * rejected attempts never contribute to the stopping coordinate.
     *
     * @return number of completed fixed output intervals.
     */
    template <class StopPredicate, class OutputCallback>
    std::size_t runUntil(
        std::size_t maximumFixedSteps,
        StopPredicate &&stopPredicate,
        OutputCallback &&outputCallback)
    {
        const auto &cfg = policy_.config();
        std::size_t completed = 0;

        for (std::size_t fixedStep = 1;
             fixedStep <= maximumFixedSteps;
             ++fixedStep)
        {
            const double targetTime =
                cfg.outputTimeOrigin +
                static_cast<double>(fixedStep) * cfg.fixedOutputDt;

            observer_.fixedTarget(fixedStep, targetTime);
            advanceTo_(targetTime);
            backend_.alignCurrentTime(targetTime);
            outputCallback(fixedStep, targetTime);
            completed = fixedStep;
            if (stopPredicate(fixedStep, targetTime))
                break;
        }
        return completed;
    }

private:
    [[noreturn]] void throwFailure_(
        const std::string &prefix,
        double targetTime,
        double attemptedDt,
        const NonlinearSolveResult &result,
        int retryCount) const
    {
        std::ostringstream message;
        message << prefix
                << " target=" << targetTime
                << " dt=" << attemptedDt
                << " retry=" << retryCount;

        if (!result.reason.empty())
            message << " reason=" << result.reason;
        else
            message << " reasonCode=" << result.reasonCode;

        throw AdaptiveTimeStepFailure(
            message.str(),
            targetTime,
            attemptedDt,
            result.reasonCode,
            retryCount);
    }

    void advanceTo_(double targetTime)
    {
        const auto &cfg = policy_.config();
        const double tolerance = cfg.effectiveTimeTolerance();
        int retryCount = 0;

        while (backend_.currentTime() < targetTime - tolerance)
        {
            const double startTime = backend_.currentTime();
            const double remaining = targetTime - startTime;
            const AdaptiveAttemptPlan plan = policy_.plan(remaining);

            backend_.setTimeStep(plan.timeStep);
            // 状态：从最近一次已接受状态建立事务快照；后续失败必须回滚到这里。
            backend_.beginAttempt();

            observer_.attempt(
                AdaptiveAttemptEvent{
                    startTime,
                    startTime + plan.timeStep,
                    plan.timeStep,
                    retryCount,
                    plan.clippedToOutputTime});

            bool nonlinearConverged = false;
            bool controlsSettled = false;
            NonlinearSolveResult lastResult;
            int stepNonlinearIterations = 0;
            long long stepLinearIterations = 0;
            double stepSolveTimeSeconds = 0.0;
            std::size_t controlResolves = 0;

            for (int controlIteration = 0;
                 controlIteration < cfg.maximumWellControlIterations;
                 ++controlIteration)
            {
                lastResult = backend_.solve();
                statistics_.recordSolve(lastResult);
                stepNonlinearIterations += lastResult.nonlinearIterations;
                stepLinearIterations += lastResult.linearIterations;
                stepSolveTimeSeconds += lastResult.wallTimeSeconds;
                nonlinearConverged = lastResult.converged;

                observer_.nonlinearSolve(
                    AdaptiveNonlinearSolveEvent{
                        startTime,
                        startTime + plan.timeStep,
                        plan.timeStep,
                        retryCount,
                        static_cast<std::size_t>(controlIteration),
                        lastResult});

                if (!nonlinearConverged)
                    break;

                // 井控制只在当前非线性系统收敛后判断；若发生切换，必须用新控制重新求解同一时间步。
                const bool switched = backend_.updateWellControls();
                if (!switched)
                {
                    controlsSettled = true;
                    break;
                }

                ++statistics_.controlSwitches;
                ++controlResolves;
                observer_.controlSwitch(
                    static_cast<std::size_t>(controlIteration + 1),
                    static_cast<std::size_t>(cfg.maximumWellControlIterations));
            }

            if (!nonlinearConverged || !controlsSettled)
            {
                // 状态：拒绝意味着解、相态和井历史都不能泄漏到下一次更小 dt 的重试。
                backend_.rejectAttempt();
                ++statistics_.rejectedSteps;

                if (!cfg.adaptive)
                {
                    throwFailure_(
                        "Fixed time step failed.",
                        targetTime,
                        plan.timeStep,
                        lastResult,
                        retryCount);
                }

                ++retryCount;
                const double nextDt = policy_.rejected(plan.timeStep);

                observer_.rejected(
                    AdaptiveRejectedEvent{
                        startTime,
                        plan.timeStep,
                        nextDt,
                        retryCount,
                        lastResult,
                        nonlinearConverged && !controlsSettled});

                if (nextDt < cfg.minimumDt)
                {
                    throwFailure_(
                        "Adaptive time step fell below minimumDt.",
                        targetTime,
                        nextDt,
                        lastResult,
                        retryCount);
                }

                if (retryCount > cfg.maximumRetries)
                {
                    throwFailure_(
                        "Maximum adaptive time-step retries exceeded.",
                        targetTime,
                        nextDt,
                        lastResult,
                        retryCount);
                }

                continue;
            }

            double acceptedTime = startTime + plan.timeStep;
            if (std::abs(acceptedTime - targetTime) <= tolerance)
                acceptedTime = targetTime;

            backend_.acceptAttempt(acceptedTime);

            statistics_.recordAcceptedStep(
                stepNonlinearIterations,
                stepLinearIterations,
                stepSolveTimeSeconds,
                plan.timeStep);

            observer_.accepted(
                AdaptiveAcceptedEvent{
                    acceptedTime,
                    plan.timeStep,
                    stepNonlinearIterations,
                    stepLinearIterations,
                    stepSolveTimeSeconds,
                    controlResolves,
                    retryCount});

            policy_.accepted(
                stepNonlinearIterations,
                plan.clippedToOutputTime,
                plan.timeStep);

            retryCount = 0;
        }
    }

    Backend &backend_;
    AdaptiveTimeStepPolicy policy_;
    Observer observer_;
    AdaptiveTimeStepStatistics statistics_;
};

} // namespace MPMC
