/**
 * @file vertical_well.hpp
 * @brief 规则竖直井穿孔构造辅助。
 */
#pragma once

#include <well/perforation.hpp>

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 在不依赖 At 或具体 Grid 的条件下生成直井穿孔。
 *
 * cellId(i,j,k) maps logical indices to the simulator's stable cell id.
 * wellIndex(i,j,k) returns the WI for that completion.
 */
template <class CellId, class CellIdFunction, class WellIndexFunction>
[[nodiscard]] std::vector<WellPerforation<CellId>>
makeVerticalPerforations(
    int i,
    int j,
    int firstLayer,
    int layerCount,
    CellIdFunction &&cellId,
    WellIndexFunction &&wellIndex)
{
    static_assert(std::is_integral_v<CellId>,
                  "Vertical perforations require integral stable cell ids.");

    if (firstLayer < 0 || layerCount <= 0)
        throw std::invalid_argument(
            "Vertical well requires a non-negative first layer and positive layer count.");

    std::vector<WellPerforation<CellId>> result;
    result.reserve(static_cast<std::size_t>(layerCount));

    for (int offset = 0; offset < layerCount; ++offset)
    {
        const int k = firstLayer + offset;
        WellPerforation<CellId> perforation{
            static_cast<CellId>(cellId(i, j, k)),
            static_cast<double>(wellIndex(i, j, k))};
        perforation.validate();
        result.push_back(perforation);
    }

    return result;
}

/**
 * @brief 根据显式给定的 WI 生成直井穿孔。
 */
template <class CellId, class CellIdFunction>
[[nodiscard]] std::vector<WellPerforation<CellId>>
makeVerticalPerforations(
    int i,
    int j,
    int firstLayer,
    const std::vector<double> &wellIndices,
    CellIdFunction &&cellId)
{
    if (wellIndices.empty())
        throw std::invalid_argument("Explicit vertical WI list cannot be empty.");

    return makeVerticalPerforations<CellId>(
        i,
        j,
        firstLayer,
        static_cast<int>(wellIndices.size()),
        std::forward<CellIdFunction>(cellId),
        [&](int, int, int k)
        {
            return wellIndices[static_cast<std::size_t>(k - firstLayer)];
        });
}

} // namespace MPMC
