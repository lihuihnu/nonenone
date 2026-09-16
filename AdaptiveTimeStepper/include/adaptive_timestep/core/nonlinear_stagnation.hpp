/**
 * @file nonlinear_stagnation.hpp
 * @brief 与 PETSc 解耦的非线性残差停滞检测器。
 */
#pragma once

#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 非线性求解残差长期无改善时的提前终止配置。
 *
 * This guard is deliberately separate from the adaptive time-step policy.  It does
 * not decide the next dt; it only prevents an already-stagnant SNES attempt from
 * spending dozens of additional Newton/KSP iterations before the stepper performs
 * the normal transactional rollback and dt cut.
 */
struct NonlinearStagnationConfig final
{
    bool enabled{false};
    int minimumIterations{8};
    int stagnantIterations{5};
    double requiredRelativeImprovement{1.0e-4};

    void validate() const
    {
        if (minimumIterations < 1)
            throw std::invalid_argument("minimumIterations must be at least 1.");
        if (stagnantIterations < 1)
            throw std::invalid_argument("stagnantIterations must be at least 1.");
        if (!std::isfinite(requiredRelativeImprovement) ||
            requiredRelativeImprovement <= 0.0 ||
            requiredRelativeImprovement >= 1.0)
        {
            throw std::invalid_argument(
                "requiredRelativeImprovement must lie in (0,1).");
        }
    }
};

/**
 * @brief 根据最近一次有意义的残差下降判断非线性平台。
 *
 * A new reference point is recorded only after the residual has dropped by at
 * least `requiredRelativeImprovement` relative to the previous reference.  Tiny
 * roundoff-level changes therefore do not keep a failed Newton solve alive.
 */
class NonlinearStagnationDetector final
{
public:
    explicit NonlinearStagnationDetector(NonlinearStagnationConfig config = {})
        : config_(config)
    {
        config_.validate();
    }

    void reset() noexcept
    {
        initialized_ = false;
        referenceResidual_ = 0.0;
        lastMeaningfulImprovementIteration_ = 0;
    }

    [[nodiscard]] const NonlinearStagnationConfig &config() const noexcept
    {
        return config_;
    }

    /** @brief 当前非线性残差已形成持续平台时返回 true。 */
    [[nodiscard]] bool update(int iteration, double residualNorm) noexcept
    {
        if (!config_.enabled || iteration < 0 ||
            !std::isfinite(residualNorm) || residualNorm < 0.0)
        {
            return false;
        }

        if (!initialized_ || iteration == 0)
        {
            initialized_ = true;
            referenceResidual_ = residualNorm;
            lastMeaningfulImprovementIteration_ = iteration;
            return false;
        }

        const double meaningfulTarget =
            referenceResidual_ * (1.0 - config_.requiredRelativeImprovement);
        if (residualNorm <= meaningfulTarget)
        {
            referenceResidual_ = residualNorm;
            lastMeaningfulImprovementIteration_ = iteration;
            return false;
        }

        if (iteration < config_.minimumIterations)
            return false;

        return iteration - lastMeaningfulImprovementIteration_ >=
               config_.stagnantIterations;
    }

private:
    NonlinearStagnationConfig config_{};
    bool initialized_{false};
    double referenceResidual_{0.0};
    int lastMeaningfulImprovementIteration_{0};
};

} // namespace MPMC
