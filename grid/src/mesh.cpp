/**
 * @file mesh.cpp
 * @brief CpGrid Mesh 高层生命周期、communicator 元数据和编号导出。
 */
#include <cpgrid/mesh.hpp>

#include "mesh_canonical_data.hpp"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>

namespace MPMC
{

Mesh::Mesh(CpCommunicator comm)
    : comm_(comm)
{
    if (comm_ == MPI_COMM_NULL)
        throw std::invalid_argument(
            "Mesh requires a valid MPI communicator.");
    initializeCommunicator_();
}

Mesh::Mesh(
    std::string dataDirectory,
    CpCommunicator comm)
    : comm_(comm),
      dataDirectory_(std::move(dataDirectory))
{
    if (comm_ == MPI_COMM_NULL)
        throw std::invalid_argument(
            "Mesh requires a valid MPI communicator.");
    initializeCommunicator_();

    if (dataDirectory_.empty())
        throw std::invalid_argument(
            "Mesh data directory is empty.");

    while (!dataDirectory_.empty() && dataDirectory_.back() == '/')
        dataDirectory_.pop_back();

    const auto data =
        detail::loadMrstCanonicalMeshData(dataDirectory_);
    if (data.cells.size() >
        static_cast<std::size_t>(
            std::numeric_limits<PetscInt>::max()))
    {
        throw std::overflow_error(
            "Mesh cell count does not fit PetscInt.");
    }
    initializeFromCanonical_(data);
}

Mesh::Mesh(
    const GrdeclGridData &grdecl,
    CpCommunicator comm)
    : comm_(comm)
{
    if (comm_ == MPI_COMM_NULL)
        throw std::invalid_argument(
            "Mesh requires a valid MPI communicator.");
    initializeCommunicator_();

    const auto data =
        detail::makeGrdeclCanonicalMeshData(grdecl);
    initializeFromCanonical_(data);
}

void Mesh::initializeCommunicator_()
{
    PetscCallMPIAbort(
        comm_,
        MPI_Comm_size(comm_, &processCount_));
    PetscCallMPIAbort(
        comm_,
        MPI_Comm_rank(comm_, &rank_));
    if (processCount_ <= 0 ||
        rank_ < 0 ||
        rank_ >= processCount_)
    {
        throw std::runtime_error(
            "Mesh communicator metadata are invalid.");
    }
}

int Mesh::processCount() const noexcept
{
    return processCount_;
}

int Mesh::rank() const noexcept
{
    return rank_;
}

void Mesh::requireTopology_() const
{
    if (!topologyInitialized_)
        throw std::logic_error(
            "Mesh topology has not been initialized.");
}

void Mesh::requirePartitioned_() const
{
    if (!partitioned_)
        throw std::logic_error(
            "Mesh has not been partitioned.");
}

void Mesh::requireCurrentIds_() const
{
    if (!currentIdsReady_)
        throw std::logic_error(
            "Mesh current ids are not available.");
}

void Mesh::writeCellIdMap(
    const std::string &filename) const
{
    requireCurrentIds_();

    if (rank() != 0)
        return;

    std::ofstream file(filename);
    if (!file)
        throw std::runtime_error(
            "Cannot open cell id map file: " + filename);

    const auto cartesianByInput = cartesianDirectory_.cartesianIdsInInputOrder();
    file << "input_index,cartesian_id,current_id,owner_rank\n";
    for (std::size_t input = 0; input < globalCellCount_; ++input)
    {
        const PetscInt currentId = currentIdFromInputIndex(input);
        file
            << input << ','
            << cartesianByInput.at(input) << ','
            << currentId << ','
            << ownerRankFromCurrentId(currentId)
            << '\n';
    }

}

} // namespace MPMC
