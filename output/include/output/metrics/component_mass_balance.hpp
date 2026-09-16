/**
 * @file component_mass_balance.hpp
 * @brief 逐组分全局质量守恒误差计算。
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace MPMC
{

/**
 * @brief 一个时刻的逐组分全局质量守恒快照。
 *
 * Sign convention follows the reservoir equations: injection is positive and
 * production is negative.  `cumulativeProduced` is stored as a positive
 * engineering magnitude, therefore the closed-boundary balance is
 *
 *   current = initial + cumulativeInjected - cumulativeProduced + error.
 */
template <std::size_t N>
struct ComponentMassBalanceSnapshot final
{
    std::array<double, N> initial{};
    std::array<double, N> current{};
    std::array<double, N> cumulativeInjected{};
    std::array<double, N> cumulativeProduced{};
    std::array<double, N> expected{};
    std::array<double, N> error{};
    std::array<double, N> relativeError{};
};

/**
 * @brief 与 PETSc 无关的 Backward-Euler 全局逐组分质量守恒账本。
 *
 * The simulator evaluates well source rates at the accepted end-of-step state.
 * That is the same time level used by the fully implicit residual, so each
 * accepted internal step contributes `rate^(n+1) * dt`.  Rejected attempts are
 * never recorded.
 */
template <std::size_t N>
class ComponentMassBalanceLedger final
{
public:
    using Array = std::array<double, N>;
    using Snapshot = ComponentMassBalanceSnapshot<N>;

    void initialize(const Array &initialInventory, double initialTimeSeconds = 0.0)
    {
        validateFinite_(initialInventory, "initial component inventory");
        if (!std::isfinite(initialTimeSeconds))
            throw std::invalid_argument("Initial mass-balance time must be finite.");

        initial_ = initialInventory;
        cumulativeInjected_.fill(0.0);
        cumulativeProduced_.fill(0.0);
        lastAcceptedTime_ = initialTimeSeconds;
        initialized_ = true;
    }

    [[nodiscard]] bool initialized() const noexcept { return initialized_; }
    [[nodiscard]] double lastAcceptedTime() const noexcept { return lastAcceptedTime_; }

    /**
     * @brief 记录一个已经接受的全隐式时间步。
     *
     * @param acceptedTimeSeconds End time of the accepted step.
     * @param injectedRate Positive per-component injection rates [kg/s].
     * @param producedRate Positive per-component production magnitudes [kg/s].
     */
    void accept(
        double acceptedTimeSeconds,
        const Array &injectedRate,
        const Array &producedRate)
    {
        requireInitialized_();
        validateFinite_(injectedRate, "component injection rate");
        validateFinite_(producedRate, "component production rate");
        if (!std::isfinite(acceptedTimeSeconds))
            throw std::invalid_argument("Accepted mass-balance time must be finite.");

        const double dt = acceptedTimeSeconds - lastAcceptedTime_;
        if (!(dt > 0.0))
            throw std::invalid_argument("Accepted mass-balance times must increase strictly.");

        for (std::size_t c = 0; c < N; ++c)
        {
            if (injectedRate[c] < -1.0e-12 || producedRate[c] < -1.0e-12)
                throw std::invalid_argument(
                    "Mass-balance injected/produced magnitudes must be non-negative.");
            cumulativeInjected_[c] += std::max(0.0, injectedRate[c]) * dt;
            cumulativeProduced_[c] += std::max(0.0, producedRate[c]) * dt;
        }
        lastAcceptedTime_ = acceptedTimeSeconds;
    }

    [[nodiscard]] Snapshot snapshot(const Array &currentInventory) const
    {
        requireInitialized_();
        validateFinite_(currentInventory, "current component inventory");

        Snapshot result;
        result.initial = initial_;
        result.current = currentInventory;
        result.cumulativeInjected = cumulativeInjected_;
        result.cumulativeProduced = cumulativeProduced_;

        for (std::size_t c = 0; c < N; ++c)
        {
            result.expected[c] =
                initial_[c] + cumulativeInjected_[c] - cumulativeProduced_[c];
            result.error[c] = currentInventory[c] - result.expected[c];

            const double scale = std::max({
                1.0,
                std::abs(initial_[c]),
                std::abs(currentInventory[c]),
                cumulativeInjected_[c] + cumulativeProduced_[c]});
            result.relativeError[c] = result.error[c] / scale;
        }
        return result;
    }

private:
    static void validateFinite_(const Array &values, const char *label)
    {
        for (double value : values)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument(std::string(label) + " must be finite.");
        }
    }

    void requireInitialized_() const
    {
        if (!initialized_)
            throw std::logic_error("Component mass-balance ledger is not initialized.");
    }

    Array initial_{};
    Array cumulativeInjected_{};
    Array cumulativeProduced_{};
    double lastAcceptedTime_{0.0};
    bool initialized_{false};
};

} // namespace MPMC
