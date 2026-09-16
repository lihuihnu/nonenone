/**
 * @file metis_graph.hpp
 * @brief 由网格邻接关系构造 METIS 分区图。
 */
#pragma once

#include <metis.h>

#include <vector>

namespace MPMC
{

class Mesh;

/**
 * @brief METIS 使用的 CSR 邻接图。
 */
struct METISCsrGraph
{
    std::vector<idx_t> offsets;
    std::vector<idx_t> neighbors;

    [[nodiscard]] idx_t vertexCount() const noexcept
    {
        return offsets.empty()
                   ? 0
                   : static_cast<idx_t>(
                         offsets.size() - 1);
    }
};

/**
 * @brief 根据 Mesh 单元邻接关系创建 METIS CSR 图。
 */
[[nodiscard]] METISCsrGraph
createMetisCsrGraph(
    const Mesh &mesh);

} // namespace MPMC
