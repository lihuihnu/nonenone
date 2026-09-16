/**
 * @file metis_partitioner.hpp
 * @brief 基于 METIS 的网格分区及 MPI 归属计算。
 */
#pragma once

#include <cpgrid/metis_graph.hpp>

namespace MPMC
{

class Mesh;

/**
 * @brief 使用 METIS 将 Mesh 单元映射到 MPI ranks。
 *
 * 当前策略固定为“一 rank 一 partition”，避免额外的 partition-to-rank 映射层。
 */
class MetisPartitioner final
{
  public:
    /**
     * @brief 执行分区并广播 owner 映射。
     *
     * `graph` 只要求在 communicator rank 0 上有效；非 root rank
     * 可以传入空 CSR，因为 METIS 仅在 root 执行。
     */
    void partition(
        Mesh &mesh,
        const METISCsrGraph &graph) const;

    /**
     * @brief 只在指定 root 上执行与 partition() 相同的 METIS 分区，不广播。
     */
    void partitionRootOnly(
        Mesh &mesh,
        const METISCsrGraph &graph,
        int root = 0) const;
};

} // namespace MPMC
