/**
 * @file mesh_local_snapshot.cpp
 * @brief CpGrid Mesh 在分区后压缩为 owned+ghost 本地拓扑快照的实现。
 */
#include <cpgrid/mesh.hpp>

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace MPMC
{

std::size_t Mesh::materializedStorageIndexFromCurrentId_(
    PetscInt currentId) const
{
    if (!localSnapshotCompacted_)
    {
        if (currentId < 0 ||
            static_cast<std::size_t>(currentId) >= inputIndexByCurrentId_.size())
        {
            throw std::out_of_range(
                "Current cell id is out of range.");
        }
        return static_cast<std::size_t>(
            inputIndexByCurrentId_[static_cast<std::size_t>(currentId)]);
    }

    if (!ownedCellIds_.empty())
    {
        const PetscInt firstOwned = ownedCellIds_.front();
        const PetscInt ownedOffset = currentId - firstOwned;
        if (ownedOffset >= 0 &&
            static_cast<std::size_t>(ownedOffset) < ownedCellIds_.size() &&
            ownedCellIds_[static_cast<std::size_t>(ownedOffset)] == currentId)
        {
            return static_cast<std::size_t>(ownedOffset);
        }
    }

    const auto ghost = std::lower_bound(
        ghostCellIds_.begin(),
        ghostCellIds_.end(),
        currentId);
    if (ghost == ghostCellIds_.end() || *ghost != currentId)
    {
        throw std::out_of_range(
            "Requested current cell is not materialized on this MPI rank.");
    }
    return ownedCellIds_.size() +
           static_cast<std::size_t>(ghost - ghostCellIds_.begin());
}

Polyhedron &Mesh::cellByInputIndex(std::size_t index)
{
    if (!localSnapshotCompacted_)
        return cells_.at(index);

    if (index >= currentIdByInputIndex_.size())
        throw std::out_of_range(
            "Requested input cell is not materialized on this MPI rank.");

    try
    {
        return cells_.at(materializedStorageIndexFromCurrentId_(
            currentIdByInputIndex_[index]));
    }
    catch (const std::out_of_range &)
    {
        throw std::out_of_range(
            "Requested input cell is not materialized on this MPI rank.");
    }
}

const Polyhedron &Mesh::cellByInputIndex(std::size_t index) const
{
    if (!localSnapshotCompacted_)
        return cells_.at(index);

    if (index >= currentIdByInputIndex_.size())
        throw std::out_of_range(
            "Requested input cell is not materialized on this MPI rank.");

    try
    {
        return cells_.at(materializedStorageIndexFromCurrentId_(
            currentIdByInputIndex_[index]));
    }
    catch (const std::out_of_range &)
    {
        throw std::out_of_range(
            "Requested input cell is not materialized on this MPI rank.");
    }
}

const Node &Mesh::node(std::size_t index) const
{
    if (!localSnapshotCompacted_)
        return nodes_.at(index);

    const PetscInt requested = static_cast<PetscInt>(index);
    const auto found = std::lower_bound(
        nodes_.begin(),
        nodes_.end(),
        requested,
        [](const Node &node, PetscInt id)
        {
            return node.id() < id;
        });
    if (found == nodes_.end() || found->id() != requested)
        throw std::out_of_range(
            "Requested node is not materialized on this MPI rank.");
    return *found;
}

void Mesh::validateMaterializedOrdering_() const
{
    if (!localSnapshotCompacted_)
        return;

    if (cells_.size() != ownedCellIds_.size() + ghostCellIds_.size())
        throw std::logic_error(
            "Mesh local snapshot cell count does not match owned+ghost ids.");

    for (std::size_t local = 0; local < ownedCellIds_.size(); ++local)
        if (cells_[local].id() != ownedCellIds_[local])
            throw std::logic_error(
                "Mesh local snapshot owned-cell ordering is inconsistent.");

    for (std::size_t ghost = 0; ghost < ghostCellIds_.size(); ++ghost)
        if (cells_[ownedCellIds_.size() + ghost].id() != ghostCellIds_[ghost])
            throw std::logic_error(
                "Mesh local snapshot ghost-cell ordering is inconsistent.");

    for (std::size_t local = 1; local < nodes_.size(); ++local)
        if (!(nodes_[local - 1].id() < nodes_[local].id()))
            throw std::logic_error(
                "Mesh local snapshot node ids are not strictly increasing.");
}

void Mesh::compactToLocalSnapshot(bool keepFullMeshOnRoot)
{
    requireCurrentIds_();

    if (localSnapshotCompacted_)
        return;

    if ((keepFullMeshOnRoot && rank() == 0) ||
        cells_.size() == ownedCellIds_.size() + ghostCellIds_.size())
    {
        return;
    }

    std::vector<PetscInt> selectedIds;
    selectedIds.reserve(ownedCellIds_.size() + ghostCellIds_.size());
    selectedIds.insert(selectedIds.end(), ownedCellIds_.begin(), ownedCellIds_.end());
    selectedIds.insert(selectedIds.end(), ghostCellIds_.begin(), ghostCellIds_.end());

    std::vector<PetscInt> selectedNodeIds;
    for (PetscInt currentId : selectedIds)
    {
        const Polyhedron &cell = cellByCurrentId(currentId);
        for (const Face &face : cell.faces())
            for (const Node *node : face.nodes())
            {
                if (node == nullptr)
                    throw std::logic_error(
                        "Mesh local snapshot encountered a null node pointer.");
                selectedNodeIds.push_back(node->id());
            }
    }

    std::sort(selectedNodeIds.begin(), selectedNodeIds.end());
    selectedNodeIds.erase(
        std::unique(selectedNodeIds.begin(), selectedNodeIds.end()),
        selectedNodeIds.end());

    std::vector<Node> localNodes;
    localNodes.reserve(selectedNodeIds.size());
    std::unordered_map<PetscInt, std::size_t> nodeLocalIndex;
    nodeLocalIndex.reserve(selectedNodeIds.size());
    for (PetscInt nodeId : selectedNodeIds)
    {
        const Node &source = nodes_.at(static_cast<std::size_t>(nodeId));
        const std::size_t local = localNodes.size();
        localNodes.emplace_back(source.x(), source.y(), source.z(), source.id());
        nodeLocalIndex.emplace(nodeId, local);
    }

    struct FaceNeighbor final
    {
        PetscInt currentId{-1};
        PetscInt inputFaceId{-1};
    };

    std::vector<Polyhedron> localCells;
    std::vector<std::vector<FaceNeighbor>> neighborByFace;
    localCells.reserve(selectedIds.size());
    neighborByFace.reserve(selectedIds.size());

    for (PetscInt currentId : selectedIds)
    {
        const Polyhedron &sourceCell = cellByCurrentId(currentId);
        std::vector<Face> faces;
        std::vector<FaceNeighbor> faceNeighbors;
        faces.reserve(sourceCell.faceCount());
        faceNeighbors.reserve(sourceCell.faceCount());

        for (const Face &sourceFace : sourceCell.faces())
        {
            std::vector<Node *> faceNodes;
            faceNodes.reserve(sourceFace.nodeCount());
            for (const Node *sourceNode : sourceFace.nodes())
            {
                const auto found = nodeLocalIndex.find(sourceNode->id());
                if (found == nodeLocalIndex.end())
                    throw std::logic_error(
                        "Mesh local snapshot lost a referenced node.");
                faceNodes.push_back(&localNodes[found->second]);
            }

            faces.emplace_back(
                std::move(faceNodes),
                sourceFace.inputPosition(),
                sourceFace.inputFaceId());

            const Polyhedron *neighbor = sourceFace.neighborCell();
            faceNeighbors.push_back({
                neighbor == nullptr ? PetscInt(-1) : neighbor->id(),
                sourceFace.inputFaceId()});
        }

        localCells.emplace_back(
            std::move(faces),
            sourceCell.inputIndex());
        Polyhedron &target = localCells.back();
        target.setId(sourceCell.id());
        target.setProcessorId(sourceCell.processorId());
        neighborByFace.push_back(std::move(faceNeighbors));
    }

    std::unordered_map<PetscInt, std::size_t> cellLocalIndex;
    cellLocalIndex.reserve(localCells.size());
    for (std::size_t local = 0; local < localCells.size(); ++local)
        cellLocalIndex.emplace(localCells[local].id(), local);

    for (std::size_t local = 0; local < localCells.size(); ++local)
    {
        Polyhedron &cell = localCells[local];
        auto &faces = cell.faces();
        for (std::size_t faceIndex = 0; faceIndex < faces.size(); ++faceIndex)
        {
            Face &face = faces[faceIndex];
            const PetscInt neighborId = neighborByFace[local][faceIndex].currentId;
            if (neighborId < 0)
                continue;

            const auto neighborFound = cellLocalIndex.find(neighborId);
            if (neighborFound == cellLocalIndex.end())
            {
                if (cell.isOwnedBy(rank()))
                    throw std::logic_error(
                        "Owned Mesh cell lost a one-ring neighbor during local compaction.");
                continue;
            }

            Polyhedron &neighborCell = localCells[neighborFound->second];
            face.setNeighborCell(&neighborCell);

            Face *opposite = nullptr;
            for (Face &candidate : neighborCell.faces())
            {
                if (candidate.inputFaceId() == face.inputFaceId())
                {
                    opposite = &candidate;
                    break;
                }
            }
            if (opposite == nullptr)
                throw std::logic_error(
                    "Mesh local snapshot cannot find the opposite face instance.");
            face.setNeighborFace(opposite);
        }
    }

    nodes_ = std::move(localNodes);
    cells_ = std::move(localCells);
    faceNeighbors_.clear();
    faceNeighbors_.shrink_to_fit();
    localSnapshotCompacted_ = true;
    validateMaterializedOrdering_();
}

} // namespace MPMC
