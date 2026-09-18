/**
 * @file policy.hpp
 * @brief 自适应时间步长的接受、拒绝与下一步推荐策略。
 */
#pragma once

#include <adaptive_timestep/core/config.hpp>
#include <adaptive_timestep/core/solve_result.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 与物理模型无关的内部时间步长推荐策略。
 *
 * Rejected attempts cut the step. Easy accepted solves grow it; difficult
 * accepted solves shrink it. The recommendation never exceeds the mandatory
 * output interval and `plan()` clips it to the next output time.
 */
class AdaptiveTimeStepPolicy final
{
public:
    explicit AdaptiveTimeStepPolicy(AdaptiveTimeStepConfig config)
        : config_(std::move(config))
    {
        config_.validate();
        recommendedDt_ = config_.effectiveMaximumDt();
    }

    /** @brief 返回已完成合法性检查的只读策略配置。 */
    [[nodiscard]] const AdaptiveTimeStepConfig &config() const noexcept
    {
        return config_;
    }

    /** @brief 给出下一次尝试步长；必要时截断到剩余输出时间。 */
    [[nodiscard]] AdaptiveAttemptPlan plan(double remainingTime) const
    {
        if (!(remainingTime > 0.0) || !std::isfinite(remainingTime))
            throw std::invalid_argument("Remaining time must be finite and positive.");

        const double tolerance = config_.effectiveTimeTolerance();
        const bool clipped = recommendedDt_ > remainingTime + tolerance;
        return AdaptiveAttemptPlan{std::min(recommendedDt_, remainingTime), clipped};
    }

    /** @brief 记录一次拒绝，并返回缩小后的重试步长。 */
    [[nodiscard]] double rejected(double actualDt)
    {
        if (!config_.adaptive)
            return recommendedDt_;

        recommendedDt_ = actualDt * config_.cutFactor;
        return recommendedDt_;
    }

    /** @brief 在内部时间步接受后更新下一步推荐值。 */
    void accepted(int nonlinearIterations, bool clipped, double actualDt)
    {
        if (!config_.adaptive)
        {
            recommendedDt_ = config_.effectiveMaximumDt();
            return;
        }

        if (clipped)
            return;

        if (nonlinearIterations <= config_.easyNonlinearIterations)
        {
            recommendedDt_ = std::min(
                config_.effectiveMaximumDt(),
                actualDt * config_.growthFactor);
        }
        else if (nonlinearIterations >= config_.difficultNonlinearIterations)
        {
            recommendedDt_ = std::max(
                config_.minimumDt,
                actualDt * config_.difficultShrinkFactor);
        }
        else
        {
            recommendedDt_ = actualDt;
        }
    }

private:
    AdaptiveTimeStepConfig config_;
    double recommendedDt_{0.0};
};

} // namespace MPMC
