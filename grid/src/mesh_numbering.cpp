/**
 * @file mesh_numbering.cpp
 * @brief Mesh current/cartesian/input 编号重排与反向查找实现。
 */
#include <cpgrid/mesh.hpp>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace MPMC
{

void Mesh::resetCurrentIds()
{
    requirePartitioned_();

    /*
     * 先统计每个 rank 的单元数，再构造 rank-major 起始偏移。
     * 第二次按稳定 storage 顺序遍历即可得到与旧 O(P*N) 双循环
     * 完全相同的 current id，同时把复杂度降为 O(N+P)。
     */
    const int ownerCount = processCount();
    std::vector<PetscInt> cellsPerOwner(
        static_cast<std::size_t>(ownerCount),
        0);

    for (const Polyhedron &cell : cells_)
    {
        const int owner = cell.processorId();
        if (owner < 0 || owner >= ownerCount)
            continue;

        ++cellsPerOwner[static_cast<std::size_t>(owner)];
    }

    rankCurrentOffsets_.assign(cellsPerOwner.size() + 1, 0);
    for (std::size_t owner = 0; owner < cellsPerOwner.size(); ++owner)
        rankCurrentOffsets_[owner + 1] =
            rankCurrentOffsets_[owner] + cellsPerOwner[owner];

    std::vector<PetscInt> nextIdByOwner(
        rankCurrentOffsets_.begin(),
        rankCurrentOffsets_.end() - 1);

    const PetscInt nextId = rankCurrentOffsets_.back();

    for (Polyhedron &cell : cells_)
    {
        const int owner = cell.processorId();
        if (owner < 0 || owner >= ownerCount)
            continue;

        auto &ownerNext =
            nextIdByOwner[static_cast<std::size_t>(owner)];
        cell.setId(ownerNext++);
    }

    if (nextId != static_cast<PetscInt>(cells_.size()))
        throw std::runtime_error(
            "resetCurrentIds did not visit every cell.");

    const PetscInt expectedLocal =
        rankCurrentOffsets_.at(static_cast<std::size_t>(rank()) + 1) -
        rankCurrentOffsets_.at(static_cast<std::size_t>(rank()));
    if (expectedLocal != localCellCount_)
        throw std::runtime_error(
            "Mesh partition count is inconsistent with current-id offsets.");

    currentIdsReady_ = true;
    rebuildCurrentIdLookup_();
    rebuildLocalPartitionView_();
}

void Mesh::rebuildCurrentIdLookup_()
{
    currentIdByInputIndex_.assign(cells_.size(), PetscInt(-1));
    inputIndexByCurrentId_.assign(cells_.size(), PetscInt(-1));

    for (std::size_t storage = 0; storage < cells_.size(); ++storage)
    {
        const PetscInt currentId = cells_[storage].id();
        if (currentId < 0 ||
            static_cast<std::size_t>(currentId) >= cells_.size())
        {
            throw std::runtime_error(
                "Cell current id is outside valid range.");
        }

        auto &inputSlot =
            inputIndexByCurrentId_[static_cast<std::size_t>(currentId)];
        if (inputSlot >= 0)
            throw std::runtime_error(
                "Duplicate current cell id.");

        currentIdByInputIndex_[storage] = currentId;
        inputSlot = static_cast<PetscInt>(storage);
    }
}

void Mesh::rebuildLocalPartitionView_()
{
    ownedCellIds_.clear();
    ghostCellIds_.clear();
    ownedCellIds_.reserve(
        static_cast<std::size_t>(localCellCount_));

    const int localRank = rank();
    for (const Polyhedron &cell : cells_)
    {
        if (cell.isOwnedBy(localRank))
            ownedCellIds_.push_back(cell.id());
    }

    if (ownedCellIds_.size() !=
        static_cast<std::size_t>(localCellCount_))
    {
        throw std::runtime_error(
            "Mesh owned-cell cache does not match localCellCount.");
    }

    /*
     * resetCurrentIds() 已保证同一 rank 内 current id 按稳定 storage 顺序连续。
     * 这里仍排序一次，使该缓存的顺序契约独立于内部遍历实现。
     */
    std::sort(
        ownedCellIds_.begin(),
        ownedCellIds_.end());

    for (PetscInt currentId : ownedCellIds_)
    {
        const Polyhedron &cell = cellByCurrentId(currentId);
        for (const Face &face : cell.faces())
        {
            const Polyhedron *neighbor = face.neighborCell();
            if (neighbor == nullptr || neighbor->isOwnedBy(localRank))
                continue;
            ghostCellIds_.push_back(neighbor->id());
        }
    }

    std::sort(
        ghostCellIds_.begin(),
        ghostCellIds_.end());
    ghostCellIds_.erase(
        std::unique(
            ghostCellIds_.begin(),
            ghostCellIds_.end()),
        ghostCellIds_.end());
}

Polyhedron &Mesh::cellByCurrentId(PetscInt currentId)
{
    requireCurrentIds_();
    if (currentId < 0 ||
        static_cast<std::size_t>(currentId) >= inputIndexByCurrentId_.size())
    {
        throw std::out_of_range(
            "Current cell id is out of range.");
    }

    return cells_.at(materializedStorageIndexFromCurrentId_(currentId));
}

const Polyhedron &Mesh::cellByCurrentId(PetscInt currentId) const
{
    requireCurrentIds_();
    if (currentId < 0 ||
        static_cast<std::size_t>(currentId) >= inputIndexByCurrentId_.size())
    {
        throw std::out_of_range(
            "Current cell id is out of range.");
    }

    return cells_.at(materializedStorageIndexFromCurrentId_(currentId));
}

Polyhedron &Mesh::cellByCartesianId(PetscInt cartesianId)
{
    return cellByInputIndex(
        static_cast<std::size_t>(cartesianDirectory_.inputIndex(cartesianId)));
}

const Polyhedron &Mesh::cellByCartesianId(PetscInt cartesianId) const
{
    return cellByInputIndex(
        static_cast<std::size_t>(cartesianDirectory_.inputIndex(cartesianId)));
}

} // namespace MPMC
