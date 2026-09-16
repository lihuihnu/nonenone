/**
 * @file format.hpp
 * @brief 输出文件中的数值、名称和列格式辅助。
 */
#pragma once

#include <output/core/types.hpp>
#include <output/metrics/component_totals.hpp>

#include <array>
#include <iomanip>
#include <ostream>
#include <string_view>
#include <vector>

namespace MPMC
{

/** @brief 写出通用逗号分隔表头 `value_0,...,value_n`。 */
inline void writeIndexedCsvHeader(std::ostream &stream, std::size_t count)
{
    for (std::size_t i = 0; i < count; ++i)
    {
        if (i != 0)
            stream << ',';
        stream << "value_" << i;
    }
    stream << '\n';
}

/** @brief 将一个采样向量写成一条 CSV 记录。 */
inline void writeCsvValues(std::ostream &stream, const std::vector<double> &values)
{
    stream << std::setprecision(16);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0)
            stream << ',';
        stream << values[i];
    }
    stream << '\n';
}

/** @brief 按 `phase_0,...,phase_n` 写出一条井相流量记录。 */
template <std::size_t NumPhases>
inline void writeWellPhaseRatesCsv(
    std::ostream &stream,
    const std::array<double, NumPhases> &phaseRates)
{
    stream << std::setprecision(16);
    for (std::size_t phase = 0; phase < NumPhases; ++phase)
    {
        if (phase != 0)
            stream << ',';
        stream << phaseRates[phase];
    }
    stream << '\n';
}

/** @brief 写出一个全局质量库存序列的标准 CSV 表头。 */
template <std::size_t NumComponents>
inline void writeMassHeader(std::ostream &stream, MassSeries series)
{
    stream << "time_seconds," << massSeriesInfo(series).totalColumn;
    for (std::size_t component = 0; component < NumComponents; ++component)
        stream << ",component_" << component << "_kg";
    stream << '\n';
}

/** @brief 按严格 CSV 格式追加一条全局质量库存记录。 */
template <std::size_t NumComponents>
inline void writeMassTotalLine(
    std::ostream &stream,
    MassSeries,
    double time,
    const ComponentTotals<NumComponents> &totals)
{
    stream << std::setprecision(16) << time << ',' << totals.total();
    for (double mass : totals.component)
        stream << ',' << mass;
    stream << '\n';
    stream.flush();
}

} // namespace MPMC
