/**
 * @file producer_composition.hpp
 * @brief 生产井逐组分累计采出、采收率和轻重选择性指标。
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace MPMC
{

template <std::size_t N>
struct ProducerCompositionSnapshot final
{
    std::array<double, N> initialInventoryKg{};
    std::array<double, N> instantaneousProducedRateKgPerS{};
    std::array<double, N> instantaneousMassFraction{};
    std::array<double, N> cumulativeProducedKg{};
    std::array<double, N> cumulativeMassFraction{};
    std::array<double, N> recoveryFraction{};
    double instantaneousLightHeavyEnrichment{
        std::numeric_limits<double>::quiet_NaN()};
    double cumulativeLightHeavyEnrichment{
        std::numeric_limits<double>::quiet_NaN()};
};

/**
 * @brief 单个生产井的 accepted-step 组分采出账本。
 *
 * Production rates are positive engineering magnitudes.  Integration follows
 * the accepted end-of-step rate, consistent with the fully implicit
 * Backward-Euler mass-balance ledger.
 */
template <std::size_t N>
class ProducerCompositionLedger final
{
public:
    using Array = std::array<double, N>;
    using Snapshot = ProducerCompositionSnapshot<N>;

    void initialize(
        const Array &initialInventoryKg,
        double initialTimeSeconds = 0.0)
    {
        validateFinite_(initialInventoryKg, "initial producer-reference inventory");
        if (!std::isfinite(initialTimeSeconds))
            throw std::invalid_argument(
                "Initial producer-composition time must be finite.");
        for (double value : initialInventoryKg)
            if (value < 0.0)
                throw std::invalid_argument(
                    "Initial producer-reference inventory must be non-negative.");

        initialInventoryKg_ = initialInventoryKg;
        cumulativeProducedKg_.fill(0.0);
        lastAcceptedTimeSeconds_ = initialTimeSeconds;
        initialized_ = true;
    }

    [[nodiscard]] bool initialized() const noexcept { return initialized_; }

    void accept(
        double acceptedTimeSeconds,
        const Array &producedRateKgPerS)
    {
        requireInitialized_();
        validateFinite_(producedRateKgPerS, "producer component rate");
        if (!std::isfinite(acceptedTimeSeconds))
            throw std::invalid_argument(
                "Accepted producer-composition time must be finite.");

        const double dt = acceptedTimeSeconds - lastAcceptedTimeSeconds_;
        if (!(dt > 0.0))
            throw std::invalid_argument(
                "Accepted producer-composition times must increase strictly.");

        for (std::size_t c = 0; c < N; ++c)
        {
            if (producedRateKgPerS[c] < -1.0e-12)
                throw std::invalid_argument(
                    "Producer component-rate magnitude must be non-negative.");
            cumulativeProducedKg_[c] +=
                std::max(0.0, producedRateKgPerS[c]) * dt;
        }
        lastAcceptedTimeSeconds_ = acceptedTimeSeconds;
    }

    [[nodiscard]] Snapshot snapshot(
        const Array &instantaneousProducedRateKgPerS,
        const std::vector<int> &lightComponents = {},
        const std::vector<int> &heavyComponents = {}) const
    {
        requireInitialized_();
        validateFinite_(
            instantaneousProducedRateKgPerS,
            "instantaneous producer component rate");

        Snapshot result;
        result.initialInventoryKg = initialInventoryKg_;
        result.instantaneousProducedRateKgPerS =
            instantaneousProducedRateKgPerS;
        result.cumulativeProducedKg = cumulativeProducedKg_;

        const double instantTotal = sum_(instantaneousProducedRateKgPerS);
        const double cumulativeTotal = sum_(cumulativeProducedKg_);
        for (std::size_t c = 0; c < N; ++c)
        {
            result.instantaneousMassFraction[c] =
                instantTotal > 0.0
                    ? instantaneousProducedRateKgPerS[c] / instantTotal
                    : std::numeric_limits<double>::quiet_NaN();
            result.cumulativeMassFraction[c] =
                cumulativeTotal > 0.0
                    ? cumulativeProducedKg_[c] / cumulativeTotal
                    : std::numeric_limits<double>::quiet_NaN();
            result.recoveryFraction[c] =
                initialInventoryKg_[c] > 0.0
                    ? cumulativeProducedKg_[c] / initialInventoryKg_[c]
                    : std::numeric_limits<double>::quiet_NaN();
        }

        result.instantaneousLightHeavyEnrichment =
            lightHeavyEnrichment_(
                result.instantaneousMassFraction,
                lightComponents,
                heavyComponents);
        result.cumulativeLightHeavyEnrichment =
            lightHeavyEnrichment_(
                result.cumulativeMassFraction,
                lightComponents,
                heavyComponents);
        return result;
    }

private:
    static void validateFinite_(const Array &values, const char *label)
    {
        for (double value : values)
            if (!std::isfinite(value))
                throw std::invalid_argument(label);
    }

    static double sum_(const Array &values)
    {
        double total = 0.0;
        for (double value : values)
            total += std::max(0.0, value);
        return total;
    }

    [[nodiscard]] double lightHeavyEnrichment_(
        const Array &producedMassFraction,
        const std::vector<int> &lightComponents,
        const std::vector<int> &heavyComponents) const
    {
        if (lightComponents.empty() || heavyComponents.empty())
            return std::numeric_limits<double>::quiet_NaN();

        const auto groupSum = [](const Array &values,
                                 const std::vector<int> &components) {
            double sum = 0.0;
            for (int component : components)
            {
                if (component < 0 ||
                    component >= static_cast<int>(N))
                    throw std::out_of_range(
                        "Producer light/heavy component index is out of range.");
                const double value =
                    values[static_cast<std::size_t>(component)];
                if (!std::isfinite(value))
                    return std::numeric_limits<double>::quiet_NaN();
                sum += value;
            }
            return sum;
        };

        const double producedLight =
            groupSum(producedMassFraction, lightComponents);
        const double producedHeavy =
            groupSum(producedMassFraction, heavyComponents);
        const double initialLight =
            groupSum(initialInventoryKg_, lightComponents);
        const double initialHeavy =
            groupSum(initialInventoryKg_, heavyComponents);

        if (!(producedHeavy > 0.0) ||
            !(initialLight > 0.0) ||
            !(initialHeavy > 0.0))
            return std::numeric_limits<double>::quiet_NaN();

        // Initial group masses and initial group mass fractions have the same
        // L/H ratio, so no separate normalization by total inventory is needed.
        return (producedLight / producedHeavy) /
               (initialLight / initialHeavy);
    }

    void requireInitialized_() const
    {
        if (!initialized_)
            throw std::logic_error(
                "ProducerCompositionLedger is not initialized.");
    }

    Array initialInventoryKg_{};
    Array cumulativeProducedKg_{};
    double lastAcceptedTimeSeconds_{0.0};
    bool initialized_{false};
};

} // namespace MPMC
