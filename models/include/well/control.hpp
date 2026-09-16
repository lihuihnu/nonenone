/**
 * @file control.hpp
 * @brief 井 BHP/率控制方程、约束和控制切换判据。
 */
#pragma once

#include <well/state.hpp>
#include <well/types.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace MPMC
{

/**
 * @brief 自动井控制使用的可选物理限制与切换滞回参数。
 *
 * 压力量单位为 Pa；流量限制保存为非负工程量，注入/采出符号由 `WellType`
 * 统一处理。
 */
struct WellControlLimits final
{
    std::optional<double> maximumBhp;
    std::optional<double> minimumBhp;
    std::optional<double> maximumWaterRate;

    double pressureTolerance{1.0e4};
    double rateRelativeTolerance{1.0e-4};
    double rateAbsoluteTolerance{1.0e-14};

    /** @brief 校验井控制限制、单位约定及非负切换滞回容差。 */
    void validate() const
    {
        const auto validatePositivePressure = [](const std::optional<double> &value,
                                                 const char *name)
        {
            if (value && (!std::isfinite(*value) || *value <= 0.0))
                throw std::invalid_argument(name);
        };

        validatePositivePressure(maximumBhp,
                                 "Maximum BHP must be positive and finite.");
        validatePositivePressure(minimumBhp,
                                 "Minimum BHP must be positive and finite.");

        if (maximumWaterRate &&
            (!std::isfinite(*maximumWaterRate) || *maximumWaterRate < 0.0))
        {
            throw std::invalid_argument(
                "Maximum water rate must be finite and non-negative.");
        }

        if (!std::isfinite(pressureTolerance) || pressureTolerance < 0.0 ||
            !std::isfinite(rateRelativeTolerance) || rateRelativeTolerance < 0.0 ||
            !std::isfinite(rateAbsoluteTolerance) || rateAbsoluteTolerance < 0.0)
        {
            throw std::invalid_argument("Well-control tolerances must be finite and non-negative.");
        }
    }
};

/** @brief 用户指定的主控制方式以及当前 Newton 方程实际采用的控制方式。 */
struct WellControlSelection final
{
    WellControl primaryControl{WellControl::Bhp};
    double primaryTarget{0.0};
    WellControl activeControl{WellControl::Bhp};
    double activeTarget{0.0};
};

/** @brief 收敛后井控制方程发生切换时记录的触发原因。 */
enum class WellControlSwitchReason
{
    None,
    MaximumBhp,
    MinimumBhp,
    PrimaryRateReached,
    PrimaryRateLimit,
    MaximumWaterRate
};

/** @brief 返回稳定的井控制切换诊断标签。 */
[[nodiscard]] constexpr std::string_view wellControlSwitchReasonName(
    WellControlSwitchReason reason) noexcept
{
    switch (reason)
    {
    case WellControlSwitchReason::None:               return "NONE";
    case WellControlSwitchReason::MaximumBhp:         return "MAX_BHP";
    case WellControlSwitchReason::MinimumBhp:         return "MIN_BHP";
    case WellControlSwitchReason::PrimaryRateReached: return "PRIMARY_RATE_REACHED";
    case WellControlSwitchReason::PrimaryRateLimit:   return "PRIMARY_RATE_LIMIT";
    case WellControlSwitchReason::MaximumWaterRate:   return "MAX_WATER_RATE";
    }
    return "UNKNOWN";
}

/**
 * @brief 一次收敛后井控制更新的可审计记录。
 *
 * `measuredValue` and `limitValue` use the physical quantity associated with
 * `reason`: pressure [Pa] for BHP limits, positive rate magnitude for rate limits.
 */
struct WellControlUpdate final
{
    bool changed{false};
    WellControl previousControl{WellControl::Bhp};
    double previousTarget{0.0};
    WellControl newControl{WellControl::Bhp};
    double newTarget{0.0};
    WellControlSwitchReason reason{WellControlSwitchReason::None};
    double measuredValue{0.0};
    double limitValue{0.0};
};

namespace detail
{

[[nodiscard]] inline bool wellRateLimitExceeded(
    double actualRate,
    double rateLimit,
    const WellControlLimits &limits) noexcept
{
    if (rateLimit <= limits.rateAbsoluteTolerance)
        return actualRate > limits.rateAbsoluteTolerance;

    return actualRate > rateLimit * (1.0 + limits.rateRelativeTolerance);
}

[[nodiscard]] inline double wellRateViolationRatio(
    double actualRate,
    double rateLimit,
    const WellControlLimits &limits) noexcept
{
    if (rateLimit <= limits.rateAbsoluteTolerance)
        return actualRate > limits.rateAbsoluteTolerance
                   ? std::numeric_limits<double>::infinity()
                   : 0.0;
    return actualRate / rateLimit;
}

} // namespace detail

/**
 * @brief 应用井控制切换规则并返回可审计记录。
 *
 * 该函数不修改储层状态，只更新 `selection`。返回的 reason/measured/limit
 * 用于日志诊断，使每次自动控制切换都能追溯其触发条件。
 */
template <class Indices>
[[nodiscard]] WellControlUpdate updateWellControl(
    WellType type,
    const WellControlLimits &limits,
    const WellState<Indices> &state,
    WellControlSelection &selection)
{
    limits.validate();

    WellControlUpdate result;
    result.previousControl = selection.activeControl;
    result.previousTarget = selection.activeTarget;
    result.newControl = selection.activeControl;
    result.newTarget = selection.activeTarget;

    if (!std::isfinite(state.bottomHolePressure))
        return result;

    const auto activate = [&](WellControl control,
                              double target,
                              WellControlSwitchReason reason,
                              double measured,
                              double limit) -> WellControlUpdate
    {
        WellControlUpdate update = result;
        if (selection.activeControl == control &&
            std::abs(selection.activeTarget - target) <= 1.0e-30)
            return update;

        update.changed = true;
        update.newControl = control;
        update.newTarget = target;
        update.reason = reason;
        update.measuredValue = measured;
        update.limitValue = limit;
        selection.activeControl = control;
        selection.activeTarget = target;
        return update;
    };

    if (type == WellType::Injector)
    {
        if (selection.activeControl != WellControl::Bhp)
        {
            if (limits.maximumBhp &&
                state.bottomHolePressure > *limits.maximumBhp + limits.pressureTolerance)
            {
                return activate(
                    WellControl::Bhp,
                    *limits.maximumBhp,
                    WellControlSwitchReason::MaximumBhp,
                    state.bottomHolePressure,
                    *limits.maximumBhp);
            }
            return result;
        }

        if (isRateControl(selection.primaryControl))
        {
            const double actual =
                state.controlledRateMagnitude(type, selection.primaryControl);
            if (detail::wellRateLimitExceeded(actual,
                                              selection.primaryTarget,
                                              limits))
            {
                return activate(
                    selection.primaryControl,
                    selection.primaryTarget,
                    WellControlSwitchReason::PrimaryRateReached,
                    actual,
                    selection.primaryTarget);
            }
        }
        return result;
    }

    // 约束：生产井最低 BHP 是最高优先级的安全限制。
    if (selection.activeControl != WellControl::Bhp &&
        limits.minimumBhp &&
        state.bottomHolePressure < *limits.minimumBhp - limits.pressureTolerance)
    {
        return activate(
            WellControl::Bhp,
            *limits.minimumBhp,
            WellControlSwitchReason::MinimumBhp,
            state.bottomHolePressure,
            *limits.minimumBhp);
    }

    const double primaryRate =
        isRateControl(selection.primaryControl)
            ? state.controlledRateMagnitude(type, selection.primaryControl)
            : 0.0;

    double waterRate = 0.0;
    if constexpr (Indices::hasWater)
        waterRate = state.controlledRateMagnitude(type, WellControl::WaterRate);

    WellControl desired = selection.activeControl;
    double desiredTarget = selection.activeTarget;
    WellControlSwitchReason desiredReason = WellControlSwitchReason::None;
    double desiredMeasured = 0.0;
    double desiredLimit = 0.0;
    double largestRatio = 1.0 + limits.rateRelativeTolerance;

    if (isRateControl(selection.primaryControl) &&
        detail::wellRateLimitExceeded(primaryRate,
                                      selection.primaryTarget,
                                      limits))
    {
        const double ratio =
            detail::wellRateViolationRatio(primaryRate,
                                           selection.primaryTarget,
                                           limits);
        if (ratio > largestRatio)
        {
            largestRatio = ratio;
            desired = selection.primaryControl;
            desiredTarget = selection.primaryTarget;
            desiredReason = WellControlSwitchReason::PrimaryRateLimit;
            desiredMeasured = primaryRate;
            desiredLimit = selection.primaryTarget;
        }
    }

    if (limits.maximumWaterRate &&
        detail::wellRateLimitExceeded(waterRate,
                                      *limits.maximumWaterRate,
                                      limits))
    {
        const double ratio =
            detail::wellRateViolationRatio(waterRate,
                                           *limits.maximumWaterRate,
                                           limits);
        if (ratio > largestRatio)
        {
            largestRatio = ratio;
            desired = WellControl::WaterRate;
            desiredTarget = *limits.maximumWaterRate;
            desiredReason = WellControlSwitchReason::MaximumWaterRate;
            desiredMeasured = waterRate;
            desiredLimit = *limits.maximumWaterRate;
        }
    }

    if (selection.activeControl == WellControl::Bhp)
    {
        if (desired != WellControl::Bhp)
            return activate(desired, desiredTarget, desiredReason, desiredMeasured, desiredLimit);
        return result;
    }

    if (desired != selection.activeControl)
        return activate(desired, desiredTarget, desiredReason, desiredMeasured, desiredLimit);

    return result;
}

} // namespace MPMC
