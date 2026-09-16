/**
 * @file metis_partitioner.cpp
 * @brief 基于 METIS 的网格分区及 MPI 归属计算的实现。
 */
#include <cpgrid/metis_partitioner.hpp>

#include <cpgrid/polyhedron.hpp>
#include <cpgrid/mesh.hpp>

#include <petscsys.h>

#include <metis.h>

#include <limits>
#include <stdexcept>
#include <vector>

namespace MPMC
{
namespace
{

void validateMetisSize(const Mesh &mesh)
{
    if (mesh.cellCount() >
        static_cast<std::size_t>(std::numeric_limits<idx_t>::max()))
        throw std::overflow_error("Mesh cell count does not fit METIS idx_t.");
}

int runMetis(
    const Mesh &mesh,
    const METISCsrGraph &graph,
    int partitionCount,
    std::vector<int> &ownerByStorage)
{
    idx_t vertexCount = static_cast<idx_t>(mesh.cellCount());
    idx_t constraintCount = 1;
    idx_t partitions = static_cast<idx_t>(partitionCount);
    idx_t edgeCut = 0;
    std::vector<idx_t> partition(mesh.cellCount(), 0);

    std::vector<idx_t> offsets = graph.offsets;
    std::vector<idx_t> neighbors = graph.neighbors;

    int status = METIS_OK;
    if (partitionCount <= 8)
    {
        status = METIS_PartGraphRecursive(
            &vertexCount,
            &constraintCount,
            offsets.data(),
            neighbors.data(),
            nullptr,
            nullptr,
            nullptr,
            &partitions,
            nullptr,
            nullptr,
            nullptr,
            &edgeCut,
            partition.data());
    }
    else
    {
        status = METIS_PartGraphKway(
            &vertexCount,
            &constraintCount,
            offsets.data(),
            neighbors.data(),
            nullptr,
            nullptr,
            nullptr,
            &partitions,
            nullptr,
            nullptr,
            nullptr,
            &edgeCut,
            partition.data());
    }

    if (status == METIS_OK)
        for (std::size_t storage = 0; storage < partition.size(); ++storage)
            ownerByStorage[storage] = static_cast<int>(partition[storage]);
    return status;
}

void applyOwners(
    Mesh &mesh,
    const std::vector<int> &owners,
    int partitionCount)
{
    if (owners.size() != mesh.cellCount())
        throw std::runtime_error("METIS owner count does not match Mesh cell count.");

    for (std::size_t storage = 0; storage < mesh.cellCount(); ++storage)
    {
        const int owner = owners[storage];
        if (owner < 0 || owner >= partitionCount)
            throw std::runtime_error("METIS produced an invalid owner rank.");
        mesh.cellByStorageIndex(storage).setProcessorId(owner);
    }
}

} // namespace

void MetisPartitioner::partition(
    Mesh &mesh,
    const METISCsrGraph &graph) const
{
    const int partitionCount = mesh.processCount();
    validateMetisSize(mesh);

    if (mesh.cellCount() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("Mesh cell count exceeds MPI_Bcast int count limit.");
    if (partitionCount <= 0)
        throw std::runtime_error("Mesh communicator has no ranks.");

    if (partitionCount == 1)
    {
        for (std::size_t storage = 0; storage < mesh.cellCount(); ++storage)
            mesh.cellByStorageIndex(storage).setProcessorId(0);
        return;
    }

    /* 保持历史 collective path 的 root CSR 校验与异常语义。 */
    if (mesh.rank() == 0 && graph.offsets.size() != mesh.cellCount() + 1)
        throw std::invalid_argument("METIS CSR offsets do not match Mesh cell count.");

    std::vector<int> ownerByStorage(mesh.cellCount(), 0);
    int metisStatus = METIS_OK;
    if (mesh.rank() == 0)
        metisStatus = runMetis(mesh, graph, partitionCount, ownerByStorage);

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Bcast(&metisStatus, 1, MPI_INT, 0, mesh.communicator()));
    if (metisStatus != METIS_OK)
        throw std::runtime_error("METIS failed to partition CpGrid.");

    const auto count = static_cast<int>(ownerByStorage.size());
    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Bcast(ownerByStorage.data(), count, MPI_INT, 0, mesh.communicator()));
    applyOwners(mesh, ownerByStorage, partitionCount);
}

void MetisPartitioner::partitionRootOnly(
    Mesh &mesh,
    const METISCsrGraph &graph,
    int root) const
{
    const int partitionCount = mesh.processCount();
    validateMetisSize(mesh);
    if (partitionCount <= 0)
        throw std::runtime_error("Mesh communicator has no ranks.");
    if (root < 0 || root >= partitionCount)
        throw std::out_of_range("METIS root rank is outside communicator.");
    if (mesh.rank() != root)
        throw std::logic_error("partitionRootOnly must be called on the selected root rank.");

    if (partitionCount == 1)
    {
        for (std::size_t storage = 0; storage < mesh.cellCount(); ++storage)
            mesh.cellByStorageIndex(storage).setProcessorId(root);
        return;
    }

    if (graph.offsets.size() != mesh.cellCount() + 1)
        throw std::invalid_argument("METIS CSR offsets do not match Mesh cell count.");

    std::vector<int> owners(mesh.cellCount(), 0);
    const int status = runMetis(mesh, graph, partitionCount, owners);
    if (status != METIS_OK)
        throw std::runtime_error("METIS failed to partition CpGrid.");
    applyOwners(mesh, owners, partitionCount);
}

} // namespace MPMC
