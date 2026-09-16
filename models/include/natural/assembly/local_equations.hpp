/**
 * @file local_equations.hpp
 * @brief 局部方程向量与自动微分 Jacobian 提取辅助。
 */
#pragma once

#include <natural/state/cell_state.hpp>

#include <array>
#include <stdexcept>

namespace MPMC
{

/** @brief 任意两热力学相之间的逸度相等残差。 */
template <class Indices, class Scalar>
[[nodiscard]] std::array<Scalar, Indices::numComponents>
fugacityEquilibriumResidual(
    const CellProperties<Indices, Scalar> &properties,
    int referencePhase,
    int otherPhase,
    double pressureScale,
    double fugacityScalingFactor)
{
    if (!(pressureScale > 0.0) || !(fugacityScalingFactor > 0.0))
        throw std::invalid_argument("Fugacity equation scaling must be positive.");
    if (referencePhase < 0 || referencePhase >= Indices::numThermodynamicPhases ||
        otherPhase < 0 || otherPhase >= Indices::numThermodynamicPhases)
        throw std::out_of_range("Thermodynamic phase index is out of range.");

    std::array<Scalar, Indices::numComponents> residual{};
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        residual[c] =
            (properties.fugacity[static_cast<std::size_t>(referencePhase)][c] -
             properties.fugacity[static_cast<std::size_t>(otherPhase)][c]) /
            pressureScale / fugacityScalingFactor;
    }
    return residual;
}

/** @brief 兼容旧油-气两相路径的快捷入口。 */
template <class Indices, class Scalar>
[[nodiscard]] std::array<Scalar, Indices::numComponents>
fugacityEquilibriumResidual(
    const CellProperties<Indices, Scalar> &properties,
    double pressureScale,
    double fugacityScalingFactor)
{
    return fugacityEquilibriumResidual<Indices>(
        properties, Indices::Phase::liquid, Indices::Phase::vapor,
        pressureScale, fugacityScalingFactor);
}

/**
 * @brief 水相 CO2 与烃相 CO2 的逸度平衡残差。
 *
 * 两相/气单相使用气相 fugacity；液单相使用液相 fugacity。
 */
template <class Indices, class Scalar>
[[nodiscard]] Scalar aqueousCO2EquilibriumResidual(
    [[maybe_unused]] const CellProperties<Indices, Scalar> &properties,
    [[maybe_unused]] HydrocarbonPhaseState phaseState,
    [[maybe_unused]] int co2Component,
    [[maybe_unused]] double pressureScale,
    [[maybe_unused]] double fugacityScalingFactor)
{
    if constexpr (!Indices::hasAqueousCO2Dissolution)
    {
        return Scalar(0.0);
    }
    else
    {
        if (co2Component < 0 || co2Component >= Indices::numComponents)
            throw std::invalid_argument("Invalid CO2 component index.");
        const std::size_t c = static_cast<std::size_t>(co2Component);
        const Scalar hydrocarbonFugacity =
            phaseState == HydrocarbonPhaseState::LiquidOnly
                ? properties.fugacity[0][c]
                : properties.fugacity[1][c];
        return (hydrocarbonFugacity - properties.aqueousCO2Fugacity) /
            pressureScale / fugacityScalingFactor;
    }
}

template <class Indices, class Scalar>
[[nodiscard]] Scalar saturationClosureResidual(
    const CellState<Indices, Scalar> &state)
{
    Scalar sum = state.liquidSaturation + state.vaporSaturation;
    if constexpr (Indices::hasWater)
        sum += state.waterSaturation;
    return sum - 1.0;
}

} // namespace MPMC
