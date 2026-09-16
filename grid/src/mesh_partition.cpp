/**
 * @file mesh_partition.cpp
 * @brief Mesh 的 METIS 分区及本地单元计数实现。
 */
#include <cpgrid/mesh.hpp>

#include <cpgrid/metis_graph.hpp>
#include <cpgrid/metis_partitioner.hpp>

#include <stdexcept>

namespace MPMC
{

void Mesh::partition()
{
    requireTopology_();

    if (partitioned_)
        throw std::logic_error(
            "Mesh::partition cannot be called twice.");

    /*
     * METIS 只在 root 上运行，因此 CSR 也只需在 root 构造。
     * 其他 rank 通过 partitioner 内部 broadcast 接收 owner 映射。
     */
    METISCsrGraph graph;
    if (processCount() > 1 && rank() == 0)
        graph = createMetisCsrGraph(*this);

    MetisPartitioner partitioner;
    partitioner.partition(*this, graph);

    localCellCount_ = 0;
    const int localRank = rank();
    for (const Polyhedron &cell : cells_)
    {
        if (cell.isOwnedBy(localRank))
            ++localCellCount_;
    }

    partitioned_ = true;
    currentIdsReady_ = false;
}

void Mesh::partitionRootOnly(int root)
{
    requireTopology_();

    if (partitioned_)
        throw std::logic_error(
            "Mesh::partitionRootOnly cannot be called twice.");
    if (root < 0 || root >= processCount())
        throw std::out_of_range(
            "Mesh root partition rank is outside communicator.");
    if (rank() != root)
        throw std::logic_error(
            "Mesh::partitionRootOnly must run on the selected root rank.");

    METISCsrGraph graph;
    if (processCount() > 1)
        graph = createMetisCsrGraph(*this);

    MetisPartitioner partitioner;
    partitioner.partitionRootOnly(*this, graph, root);

    localCellCount_ = 0;
    for (const Polyhedron &cell : cells_)
        if (cell.isOwnedBy(root))
            ++localCellCount_;

    partitioned_ = true;
    currentIdsReady_ = false;
}

} // namespace MPMC
