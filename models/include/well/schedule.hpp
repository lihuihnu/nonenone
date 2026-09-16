/**
 * @file schedule.hpp
 * @brief 井启停和阶段性控制的时间调度。
 */
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 包含端点的简单井启用时间区间。
 *
 * Preserves the old natural semantics:
 * `openTime <= time <= closeTime`.
 */
struct WellSchedule final
{
    bool enabled{true};
    double openTime{0.0};
    double closeTime{std::numeric_limits<double>::max()};

    void validate() const
    {
        if (!std::isfinite(openTime) || !std::isfinite(closeTime))
            throw std::invalid_argument("Well schedule times must be finite.");
        if (closeTime < openTime)
            throw std::invalid_argument("Well close time cannot precede open time.");
    }

    [[nodiscard]] bool isActive(double time) const noexcept
    {
        return enabled && time >= openTime && time <= closeTime;
    }
};

} // namespace MPMC
