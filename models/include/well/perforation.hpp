/**
 * @file perforation.hpp
 * @brief 单个井穿孔的网格位置、连接和局部参数。
 */
#pragma once

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace MPMC
{

/**
 * @brief 一条井与网格单元之间的穿孔连接。
 *
 * `currentCellId` is a stable global/current cell id understood by the chosen
 * grid adapter. `wellIndex` is the Peaceman connection index WI used in
 * `q = -WI * lambda * Delta p`.
 */
template <class CellId>
struct WellPerforation final
{
    static_assert(std::is_integral_v<CellId>,
                  "WellPerforation uses stable integral grid cell ids.");

    CellId currentCellId{-1};
    double wellIndex{0.0};

    /** @brief 拒绝负单元 id、负 WI、NaN 和无穷值。 */
    void validate() const
    {
        if (currentCellId < 0)
            throw std::invalid_argument("Well perforation cell id must be non-negative.");
        if (!(wellIndex >= 0.0) || !std::isfinite(wellIndex))
            throw std::invalid_argument("Well index must be finite and non-negative.");
    }
};

} // namespace MPMC
