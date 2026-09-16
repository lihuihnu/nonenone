/**
 * @file component_totals.hpp
 * @brief 各相与全域逐组分摩尔/质量总量统计。
 */
#pragma once

#include <array>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace MPMC
{

template <std::size_t NumComponents>
struct ComponentTotals final
{
    std::array<double, NumComponents> component{};

    [[nodiscard]] double total() const noexcept
    {
        return std::accumulate(
            component.begin(),
            component.end(),
            0.0);
    }
};

/**
 * @brief 对以 cell-major 形式存储的组分字段求和：
 *        [cell0 c0..cN-1][cell1 c0..cN-1]...
 */
template <std::size_t NumComponents>
[[nodiscard]] ComponentTotals<NumComponents>
sumCellMajorComponents(
    const std::vector<double> &values)
{
    static_assert(
        NumComponents > 0,
        "Component total requires at least one component.");

    if (values.size() % NumComponents != 0)
    {
        throw std::invalid_argument(
            "Cell-major field size is not divisible by component count.");
    }

    ComponentTotals<NumComponents> result;

    for (std::size_t offset = 0;
         offset < values.size();
         offset += NumComponents)
    {
        for (std::size_t component = 0;
             component < NumComponents;
             ++component)
        {
            result.component[component] +=
                values[offset + component];
        }
    }

    return result;
}

} // namespace MPMC
