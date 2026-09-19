/**
 * @file config.hpp
 * @brief 模块配置参数、默认值及输入合法性约束。
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 自适应隐式时间步策略参数。
 *
 * All time quantities use the simulator base unit (seconds in current cases).
 * `fixedOutputDt` is the distance between mandatory output times; adaptive
 * internal steps are clipped so that accepted states land exactly on them.
 */
struct AdaptiveTimeStepConfig final
{
    double fixedOutputDt{1.0};          ///< Mandatory output-time spacing.
    // Optional internal-step ceiling. Zero means "same as fixedOutputDt".
    // This separates reporting cadence from the nonlinear stability limit.
    double maximumDt{0.0};
    double minimumDt{1.0 / 1024.0};    ///< Smallest retry step accepted by the policy.
    double cutFactor{0.5};              ///< `dt_new = cutFactor * dt` after rejection.
    double growthFactor{2.0};           ///< Growth after an easy accepted solve.
    double difficultShrinkFactor{0.8}; ///< Shrink after a difficult accepted solve.

    int easyNonlinearIterations{10};       ///< SNES iterations at/below which a step is easy.
    int difficultNonlinearIterations{20};  ///< SNES iterations at/above which a step is difficult.
    int maximumRetries{10};                ///< Maximum rejected attempts at one physical time.
    int maximumWellControlIterations{6};   ///< Maximum solve/control-switch cycles per attempt.

    bool adaptive{true};
    double timeTolerance{0.0};  ///< Zero selects a scale-aware default tolerance.
    double outputTimeOrigin{0.0}; ///< Reference origin of the output-time lattice.

    /** @brief 返回实际内部时间步上限；0 表示沿用固定输出间隔。 */
    [[nodiscard]] double effectiveMaximumDt() const noexcept
    {
        return maximumDt > 0.0
            ? std::min(maximumDt, fixedOutputDt)
            : fixedOutputDt;
    }

    /** @brief 比较物理时间和输出时间时使用的数值容差。 */
    [[nodiscard]] double effectiveTimeTolerance() const noexcept
    {
        if (timeTolerance > 0.0)
            return timeTolerance;
        return std::max(1.0e-10 * fixedOutputDt, 1.0e-8);
    }

    /** @brief 检查缩放因子、阈值和时间尺度是否自洽。 */
    void validate() const
    {
        const auto finitePositive = [](double value) {
            return std::isfinite(value) && value > 0.0;
        };

        if (!finitePositive(fixedOutputDt))
            throw std::invalid_argument("fixedOutputDt must be finite and positive.");

        if (maximumDt < 0.0 || !std::isfinite(maximumDt))
            throw std::invalid_argument("maximumDt must be finite and non-negative.");

        if (!finitePositive(minimumDt) || minimumDt > effectiveMaximumDt())
            throw std::invalid_argument("minimumDt must be finite, positive, and no larger than the effective maximumDt.");

        if (!std::isfinite(cutFactor) || cutFactor <= 0.0 || cutFactor >= 1.0)
            throw std::invalid_argument("cutFactor must lie in (0,1).");

        if (!std::isfinite(growthFactor) || growthFactor < 1.0)
            throw std::invalid_argument("growthFactor must be finite and no smaller than 1.");

        if (!std::isfinite(difficultShrinkFactor) ||
            difficultShrinkFactor <= 0.0 ||
            difficultShrinkFactor > 1.0)
        {
            throw std::invalid_argument("difficultShrinkFactor must lie in (0,1].");
        }

        if (easyNonlinearIterations < 0 ||
            difficultNonlinearIterations < easyNonlinearIterations)
        {
            throw std::invalid_argument("Invalid nonlinear-iteration thresholds.");
        }

        if (maximumRetries < 1)
            throw std::invalid_argument("maximumRetries must be at least 1.");

        if (maximumWellControlIterations < 1)
            throw std::invalid_argument("maximumWellControlIterations must be at least 1.");

        if (timeTolerance < 0.0 || !std::isfinite(timeTolerance))
            throw std::invalid_argument("timeTolerance must be finite and non-negative.");

        if (!std::isfinite(outputTimeOrigin))
            throw std::invalid_argument("outputTimeOrigin must be finite.");
    }
};

} // namespace MPMC
