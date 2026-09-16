/**
 * @file property_conversion.hpp
 * @brief 运行时物性数据与 PETSc/标量存储之间的转换辅助。
 */
#pragma once

#include <common/math.hpp>
#include <natural/state/cell_state.hpp>

#include <array>
#include <cstddef>

namespace MPMC
{

/**
 * @brief 将当前 AD CellProperties 复制成上一时间层使用的纯 double 快照。
 *
 * 上一时间层蓄积项必须是 Newton 常量，因此这里只保留函数值，不保留任何导数。
 */
template <class Indices, class Scalar>
[[nodiscard]] CellProperties<Indices, double>
scalarizeCellProperties(
    const CellProperties<Indices, Scalar> &input)
{
    CellProperties<Indices, double> result;

    for (int phase = 0;
         phase < Indices::numPhases;
         ++phase)
    {
        const auto p =
            static_cast<std::size_t>(phase);

        result.density[p] =
            scalarValue(input.density[p]);
        result.viscosity[p] =
            scalarValue(input.viscosity[p]);
        result.mobility[p] =
            scalarValue(input.mobility[p]);
        result.saturation[p] =
            scalarValue(input.saturation[p]);

        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            const auto c =
                static_cast<std::size_t>(component);

            result.massFraction[p][c] =
                scalarValue(
                    input.massFraction[p][c]);
        }
    }

    for (int hydrocarbonPhase = 0;
         hydrocarbonPhase < 2;
         ++hydrocarbonPhase)
    {
        const auto p =
            static_cast<std::size_t>(
                hydrocarbonPhase);

        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            const auto c =
                static_cast<std::size_t>(component);

            result.fugacity[p][c] =
                scalarValue(
                    input.fugacity[p][c]);
        }
    }

    result.porosity =
        scalarValue(input.porosity);

    result.aqueousCO2MassFraction =
        scalarValue(
            input.aqueousCO2MassFraction);

    result.aqueousCO2Fugacity =
        scalarValue(
            input.aqueousCO2Fugacity);

    result.trappedGasSaturation =
        scalarValue(
            input.trappedGasSaturation);

    result.freeGasSaturation =
        scalarValue(
            input.freeGasSaturation);

    for (int component = 0;
         component < Indices::numComponents;
         ++component)
    {
        const auto c =
            static_cast<std::size_t>(component);

        result.adsorbedVolume[c] =
            scalarValue(
                input.adsorbedVolume[c]);
    }

    return result;
}

} // namespace MPMC
