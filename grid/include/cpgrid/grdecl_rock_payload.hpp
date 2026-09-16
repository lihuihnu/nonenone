/**
 * @file grdecl_rock_payload.hpp
 * @brief 从完整 GRDECL 数据中提取可独立保留的轻量岩石属性载荷。
 */
#pragma once

#include <cpgrid/grdecl.hpp>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief GRDECL Mesh 已构造后仍需保留的 active-cell 岩石属性。
 *
 * 该结构不保存 ACTNUM、corner geometry 或邻接信息，允许 root 在 topology
 * 构造完成后尽早释放 `GrdeclGridData::activeCells` 等大型临时数据。
 */
struct CpGridRootRockData final
{
    std::size_t cellCount{0};
    std::vector<double> porosity;
    std::array<std::vector<double>, 3> permeability;
};

/**
 * @brief 移动 PORO/PERM 到轻量载荷；不复制 expanded corner geometry。
 */
[[nodiscard]] inline CpGridRootRockData takeRootRockData(
    GrdeclGridData &grdecl)
{
    CpGridRootRockData result;
    result.cellCount = grdecl.activeCells.size();
    result.porosity = std::move(grdecl.porosity);
    result.permeability = {{
        std::move(grdecl.permeabilityX),
        std::move(grdecl.permeabilityY),
        std::move(grdecl.permeabilityZ)}};
    return result;
}

} // namespace MPMC
