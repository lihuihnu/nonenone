/**
 * @file accumulation.hpp
 * @brief 组分守恒方程的孔隙体积蓄积项计算。
 */
#pragma once

#include <natural/state/cell_state.hpp>

#include <array>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 单元流体质量蓄积量，单位为 kg/m^3 bulk volume。
 *
 * 这些量直接进入组分质量守恒方程的时间项；乘以单元 bulk volume
 * 后即得到该单元内对应守恒量的质量 [kg]。
 */
template <class Indices, class Scalar>
struct AccumulationResult
{
    std::array<Scalar, Indices::numComponents> componentMass{};
    Scalar waterMass{0.0};
};

/** @brief 直接计算单个守恒组分的流体蓄积，避免热循环构造整份临时数组。 */
template <class Indices, class Scalar>
[[nodiscard]] Scalar computeFluidComponentAccumulation(
    const CellProperties<Indices, Scalar> &properties,
    int component,
    int dissolvedCO2Component = -1)
{
    if (component < 0 || component >= Indices::numComponents)
        throw std::out_of_range("Fluid accumulation component index is out of range.");

    const std::size_t c = static_cast<std::size_t>(component);
    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        Scalar mass = 0.0;
        for (int phase = 0; phase < Indices::numPhases; ++phase)
        {
            const std::size_t p = static_cast<std::size_t>(phase);
            mass += properties.density[p] * properties.saturation[p] *
                    properties.massFraction[p][c];
        }
        return properties.porosity * mass;
    }
    else
    {
        const std::size_t liquid = static_cast<std::size_t>(Indices::Phase::liquid);
        const std::size_t vapor = static_cast<std::size_t>(Indices::Phase::vapor);
        Scalar result = properties.porosity *
            (properties.density[liquid] * properties.saturation[liquid] *
                 properties.massFraction[liquid][c] +
             properties.density[vapor] * properties.saturation[vapor] *
                 properties.massFraction[vapor][c]);

        if constexpr (Indices::hasIndependentWaterConservation &&
                      Indices::hasAqueousCO2Dissolution)
        {
            if (dissolvedCO2Component < 0 ||
                dissolvedCO2Component >= Indices::numComponents)
                throw std::invalid_argument(
                    "A valid CO2 component index is required for aqueous dissolution.");
            if (component == dissolvedCO2Component)
            {
                const std::size_t water =
                    static_cast<std::size_t>(Indices::Phase::water);
                result += properties.porosity * properties.density[water] *
                          properties.saturation[water] *
                          properties.aqueousCO2MassFraction;
            }
        }
        return result;
    }
}

/** @brief 直接计算独立水守恒蓄积；未启用独立水方程时返回零。 */
template <class Indices, class Scalar>
[[nodiscard]] Scalar computeFluidWaterAccumulation(
    const CellProperties<Indices, Scalar> &properties)
{
    if constexpr (!Indices::hasIndependentWaterConservation)
    {
        return Scalar(0.0);
    }
    else
    {
        const std::size_t water = static_cast<std::size_t>(Indices::Phase::water);
        const Scalar totalAqueousMass =
            properties.porosity * properties.density[water] *
            properties.saturation[water];
        if constexpr (Indices::hasAqueousCO2Dissolution)
            return totalAqueousMass * (1.0 - properties.aqueousCO2MassFraction);
        else
            return totalAqueousMass;
    }
}

/**
 * @brief 计算各守恒方程的流体质量蓄积项。
 *
 * 对烃类组分 i：
 * `M_i/V_b = phi (rho_l S_l w_li + rho_g S_g w_gi)`。
 * 其中 `phi` 为孔隙度，`rho` 为储层条件质量密度，`S` 为相饱和度，
 * `w` 为相内质量分数。
 *
 * 显式水相存在时：
 * `M_w/V_b = phi rho_w S_w`。启用水相 CO2 溶解后再按质量分数
 * `X_CO2,w` 拆分为
 * `M_CO2,w/V_b = phi rho_w S_w X_CO2,w` 和
 * `M_H2O/V_b = phi rho_w S_w (1-X_CO2,w)`。
 */
template <class Indices, class Scalar>
[[nodiscard]] AccumulationResult<Indices, Scalar>
computeFluidAccumulation(
    const CellProperties<Indices, Scalar> &properties,
    int dissolvedCO2Component = -1)
{
    AccumulationResult<Indices, Scalar> result;
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        result.componentMass[static_cast<std::size_t>(component)] =
            computeFluidComponentAccumulation<Indices>(
                properties, component, dissolvedCO2Component);
    }
    if constexpr (Indices::hasIndependentWaterConservation)
        result.waterMass = computeFluidWaterAccumulation<Indices>(properties);
    return result;
}

/**
 * @brief 竞争吸附造成的组分质量蓄积，单位为 kg/m^3 bulk volume。
 *
 * 第 i 个组分的吸附质量项采用
 * `M_ads,i/V_b = (1-phi) rho_r rho_i,sc theta_i`，
 * 其中 `rho_r` 为岩石骨架密度，`rho_i,sc` 为组分标准状况气体密度，
 * `theta_i` 为竞争 Langmuir 模型给出的吸附标准体积/岩石体积。
 */
template <class Indices, class Scalar>
[[nodiscard]] Scalar computeAdsorbedComponentAccumulation(
    const Scalar &porosity,
    double rockDensity,
    const std::array<double, Indices::numComponents> &standardGasDensity,
    const std::array<Scalar, Indices::numComponents> &adsorbedVolume,
    int component)
{
    if (component < 0 || component >= Indices::numComponents)
        throw std::out_of_range("Adsorbed accumulation component index is out of range.");
    const std::size_t c = static_cast<std::size_t>(component);
    return (1.0 - porosity) * rockDensity *
           standardGasDensity[c] * adsorbedVolume[c];
}

template <class Indices, class Scalar>
[[nodiscard]] std::array<Scalar, Indices::numComponents>
computeAdsorbedAccumulation(
    const Scalar &porosity,
    double rockDensity,
    const std::array<double, Indices::numComponents> &standardGasDensity,
    const std::array<Scalar, Indices::numComponents> &adsorbedVolume)
{
    std::array<Scalar, Indices::numComponents> result{};
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        result[static_cast<std::size_t>(component)] =
            computeAdsorbedComponentAccumulation<Indices>(
                porosity, rockDensity, standardGasDensity, adsorbedVolume, component);
    }
    return result;
}

} // namespace MPMC
