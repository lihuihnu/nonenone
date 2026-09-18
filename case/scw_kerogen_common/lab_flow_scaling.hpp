#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

/**
 * @file lab_flow_scaling.hpp
 * @brief 实验 slab 的 PVI、Darcy 流量/压降与井控尺度关系。
 *
 * 本文件不提供任何默认实验尺寸、渗透率、黏度或井压。
 * 所有输入必须来自装置、试件和已通过的物性验证。
 */
namespace ScwKerogenLab
{

struct FlowScale final
{
    double effectivePoreVolumeM3{0.0};
    double effectiveFlowAreaM2{0.0};
    double flowLengthM{0.0};
    double effectivePorosity{0.0};
    double absolutePermeabilityM2{0.0};
    double referenceViscosityPaS{0.0};

    void validate() const
    {
        const auto positiveFinite = [](double value) {
            return std::isfinite(value) && value > 0.0;
        };
        if (!positiveFinite(effectivePoreVolumeM3) ||
            !positiveFinite(effectiveFlowAreaM2) ||
            !positiveFinite(flowLengthM) ||
            !positiveFinite(effectivePorosity) ||
            effectivePorosity > 1.0 ||
            !positiveFinite(absolutePermeabilityM2) ||
            !positiveFinite(referenceViscosityPaS))
        {
            throw std::invalid_argument(
                "Laboratory flow-scale inputs must be finite, positive, and physical.");
        }
    }

    [[nodiscard]] double injectionRateFromPviRate(double pviPerSecond) const
    {
        validate();
        if (!(std::isfinite(pviPerSecond) && pviPerSecond > 0.0))
            throw std::invalid_argument("PVI rate must be finite and positive.");
        return pviPerSecond * effectivePoreVolumeM3;
    }

    [[nodiscard]] double nominalResidenceTime(double reservoirRateM3PerS) const
    {
        validate();
        if (!(std::isfinite(reservoirRateM3PerS) && reservoirRateM3PerS > 0.0))
            throw std::invalid_argument("Reservoir injection rate must be finite and positive.");
        return effectivePoreVolumeM3 / reservoirRateM3PerS;
    }

    [[nodiscard]] double darcyVelocity(double reservoirRateM3PerS) const
    {
        validate();
        if (!(std::isfinite(reservoirRateM3PerS) && reservoirRateM3PerS > 0.0))
            throw std::invalid_argument("Reservoir injection rate must be finite and positive.");
        return reservoirRateM3PerS / effectiveFlowAreaM2;
    }

    [[nodiscard]] double poreVelocity(double reservoirRateM3PerS) const
    {
        return darcyVelocity(reservoirRateM3PerS) / effectivePorosity;
    }

    [[nodiscard]] double singlePhaseDarcyPressureDrop(
        double reservoirRateM3PerS) const
    {
        validate();
        return referenceViscosityPaS * flowLengthM *
            darcyVelocity(reservoirRateM3PerS) /
            absolutePermeabilityM2;
    }
};

[[nodiscard]] inline double actualPvi(
    double cumulativeInjectedReservoirVolumeM3,
    double effectivePoreVolumeM3)
{
    if (!std::isfinite(cumulativeInjectedReservoirVolumeM3) ||
        cumulativeInjectedReservoirVolumeM3 < 0.0 ||
        !(std::isfinite(effectivePoreVolumeM3) && effectivePoreVolumeM3 > 0.0))
    {
        throw std::invalid_argument(
            "PVI requires non-negative cumulative injection and positive finite pore volume.");
    }
    return cumulativeInjectedReservoirVolumeM3 / effectivePoreVolumeM3;
}

[[nodiscard]] inline double producerBhpFromCenteredPressure(
    double pressureCenterPa,
    double designPressureDropPa)
{
    if (!(std::isfinite(pressureCenterPa) && pressureCenterPa > 0.0) ||
        !(std::isfinite(designPressureDropPa) && designPressureDropPa >= 0.0))
        throw std::invalid_argument("Pressure-center design inputs are invalid.");
    return pressureCenterPa - 0.5 * designPressureDropPa;
}

[[nodiscard]] inline double expectedInjectorPressureFromCenteredPressure(
    double pressureCenterPa,
    double designPressureDropPa)
{
    if (!(std::isfinite(pressureCenterPa) && pressureCenterPa > 0.0) ||
        !(std::isfinite(designPressureDropPa) && designPressureDropPa >= 0.0))
        throw std::invalid_argument("Pressure-center design inputs are invalid.");
    return pressureCenterPa + 0.5 * designPressureDropPa;
}

[[nodiscard]] inline double injectorMaximumBhp(
    double apparatusSafeLimitPa,
    double validatedPvtUpperBoundPa,
    double experimentDesignCapPa)
{
    const auto positiveFinite = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (!positiveFinite(apparatusSafeLimitPa) ||
        !positiveFinite(validatedPvtUpperBoundPa) ||
        !positiveFinite(experimentDesignCapPa))
    {
        throw std::invalid_argument("Injector BHP limits must be finite and positive.");
    }
    return std::min({
        apparatusSafeLimitPa,
        validatedPvtUpperBoundPa,
        experimentDesignCapPa});
}

} // namespace ScwKerogenLab
