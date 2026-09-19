/**
 * @file state.hpp
 * @brief 井 BHP、相流量、组分流量及控制状态。
 */
#pragma once

#include <well/types.hpp>

#include <array>
#include <cstddef>
#include <limits>

namespace MPMC
{

/**
 * @brief 由当前收敛储层状态计算得到的井状态量。
 *
 * 各相流量均采用模拟器符号约定（注入为正、采出为负）。`surfacePhaseRate`
 * 用于流量控制；`reservoirPhaseRate`、`phaseMassRate` 和
 * `flowWeightedDensity` 主要用于诊断和输出。
 */
template <class Indices>
struct WellState final
{
    double bottomHolePressure{std::numeric_limits<double>::quiet_NaN()};
    std::array<double, Indices::numPhases> surfacePhaseRate{};
    std::array<double, Indices::numPhases> reservoirPhaseRate{};
    std::array<double, Indices::numPhases> phaseMassRate{};
    std::array<std::array<double, Indices::numComponents>, Indices::numPhases>
        phaseComponentMassRate{};
    std::array<double, Indices::numPhases> flowWeightedDensity{};
    /**
     * Signed conserved-component mass rates [kg/s].
     *
     * The sign follows the simulator convention: injection is positive and
     * production is negative.  Unlike a phase-composition sample at the well
     * cell, these rates are accumulated from the exact perforation source
     * kernel and therefore remain meaningful when several phases flow at once.
     */
    std::array<double, Indices::numComponents> componentMassRate{};

    /** @brief 清零流量/密度累积量，但保留当前 BHP。 */
    void clearRates() noexcept
    {
        surfacePhaseRate.fill(0.0);
        reservoirPhaseRate.fill(0.0);
        phaseMassRate.fill(0.0);
        for (auto &phaseRate : phaseComponentMassRate)
            phaseRate.fill(0.0);
        flowWeightedDensity.fill(0.0);
        componentMassRate.fill(0.0);
    }

    /** @brief 返回指定流量控制对应的带符号地面流量。 */
    [[nodiscard]] double signedControlledRate(WellControl control) const
    {
        switch (control)
        {
        case WellControl::TotalRate:
        {
            double total = 0.0;
            for (double rate : surfacePhaseRate) total += rate;
            return total;
        }
        case WellControl::ReservoirTotalRate:
        {
            double total = 0.0;
            for (double rate : reservoirPhaseRate) total += rate;
            return total;
        }
        case WellControl::OilRate:
            return surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::liquid)];
        case WellControl::GasRate:
            return surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::vapor)];
        case WellControl::WaterRate:
            if constexpr (Indices::hasWater)
                return surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::water)];
            return 0.0;
        case WellControl::Bhp:
            return 0.0;
        }
        return 0.0;
    }

    /** @brief 返回指定控制流量的非负工程量。 */
    [[nodiscard]] double controlledRateMagnitude(WellType type, WellControl control) const
    {
        return rateMagnitude(type, signedControlledRate(control));
    }
};

} // namespace MPMC
