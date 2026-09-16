/**
 * @file dof_map.cpp
 * @brief 分区网格单元与 PETSc 全局/局部自由度编号映射的实现。
 */
#include <cpgrid/dof_map.hpp>

#include <petscsys.h>

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace MPMC
{

DofMap::DofMap(
    Mesh &mesh,
    PetscInt dofPerCell)
    : mesh_(mesh),
      dofPerCell_(dofPerCell)
{
    if (!mesh_.isReadyForDofs())
    {
        throw std::invalid_argument(
            "DofMap requires a Mesh prepared for DOF creation.");
    }

    if (dofPerCell_ <= 0)
    {
        throw std::invalid_argument(
            "DofMap requires dofPerCell > 0.");
    }

    buildOwnership_();
    buildGhostMap_();
}

void DofMap::buildOwnership_()
{
    const PetscInt globalCells =
        static_cast<PetscInt>(
            mesh_.cellCount());

    if (globalCells >
        std::numeric_limits<PetscInt>::max() /
            dofPerCell_)
    {
        throw std::overflow_error(
            "Global CpGrid DOF count overflows PetscInt.");
    }

    globalDofCount_ =
        globalCells *
        dofPerCell_;

    ownedCellCount_ =
        mesh_.localCellCount();

    firstOwnedCellId_ = mesh_.firstOwnedCurrentId();

    /*
     * Mesh::resetCurrentIds() 必须保证每个 rank 的 current ids 连续。
     * G8 开始直接消费 Mesh 已缓存的 owned-id 视图，避免每注册一种
     * DOF layout 都重新扫描完整 replicated cell storage。
     */
    const auto &ownedCellIds = mesh_.ownedCellIds();
    if (ownedCellIds.size() !=
        static_cast<std::size_t>(ownedCellCount_))
    {
        throw std::runtime_error(
            "DofMap local cell count is inconsistent with Mesh.");
    }

    for (PetscInt localOrdinal = 0;
         localOrdinal < ownedCellCount_;
         ++localOrdinal)
    {
        if (ownedCellIds[static_cast<std::size_t>(localOrdinal)] !=
            firstOwnedCellId_ + localOrdinal)
        {
            throw std::runtime_error(
                "Mesh current cell ids are not contiguous inside a rank.");
        }
    }
}

void DofMap::buildGhostMap_()
{
    /*
     * ghost 拓扑属于 Mesh 分区本身，与具体 DOF 布局无关。DofMap 只引用
     * Mesh 的稳定 ghost 目录，不再为每个 layout 复制一份 O(Nghost) 数组。
     */
    const auto &ghostCellIds = mesh_.ghostCellIds();
    if (!std::is_sorted(ghostCellIds.begin(), ghostCellIds.end()))
        throw std::logic_error("Mesh ghost cell ids must be sorted.");
}

void DofMap::validateComponent_(
    PetscInt component) const
{
    if (component < 0 ||
        component >= dofPerCell_)
    {
        throw std::out_of_range(
            "DOF component is outside [0, dofPerCell).");
    }
}

PetscInt DofMap::globalIndex(
    const Polyhedron &cell,
    PetscInt component) const
{
    return globalIndex(cell.id(), component);
}

PetscInt DofMap::globalIndex(
    PetscInt currentCellId,
    PetscInt component) const
{
    validateComponent_(component);

    const PetscInt globalCellCount =
        globalDofCount_ / dofPerCell_;
    if (currentCellId < 0 || currentCellId >= globalCellCount)
        throw std::out_of_range(
            "Computed global DOF is outside global range.");

    return currentCellId * dofPerCell_ + component;
}

PetscInt DofMap::localBlockIndex(
    const Polyhedron &cell) const
{
    return localBlockIndex(cell.id());
}

PetscInt DofMap::localBlockIndex(
    PetscInt currentCellId) const
{
    if (isOwnedCellId(currentCellId))
        return currentCellId - firstOwnedCellId_;

    const auto found = std::lower_bound(
        ghostCellIds().begin(),
        ghostCellIds().end(),
        currentCellId);
    if (found == ghostCellIds().end() || *found != currentCellId)
    {
        throw std::out_of_range(
            "Remote cell is not present in the current ghost stencil.");
    }
    return ownedCellCount_ +
           static_cast<PetscInt>(found - ghostCellIds().begin());
}

PetscInt DofMap::localIndex(
    const Polyhedron &cell,
    PetscInt component) const
{
    return localIndex(cell.id(), component);
}

PetscInt DofMap::localIndex(
    PetscInt currentCellId,
    PetscInt component) const
{
    validateComponent_(component);
    return localBlockIndex(currentCellId) * dofPerCell_ + component;
}

std::vector<PetscInt>
DofMap::localToGlobalBlockIndices() const
{
    std::vector<PetscInt> blocks;
    blocks.reserve(
        static_cast<std::size_t>(
            ownedCellCount_) +
        ghostCellIds().size());

    for (PetscInt local = 0;
         local < ownedCellCount_;
         ++local)
    {
        blocks.push_back(
            firstOwnedCellId_ +
            local);
    }

    blocks.insert(
        blocks.end(),
        ghostCellIds().begin(),
        ghostCellIds().end());

    return blocks;
}

std::vector<PetscInt>
DofMap::localToGlobalIndices() const
{
    const auto &ghostIds = ghostCellIds();
    const std::size_t localCellCount =
        static_cast<std::size_t>(ownedCellCount_) +
        ghostIds.size();

    std::vector<PetscInt> indices;
    indices.reserve(
        localCellCount *
        static_cast<std::size_t>(dofPerCell_));

    const auto appendCell =
        [this, &indices](PetscInt block)
        {
            const PetscInt first = block * dofPerCell_;
            for (PetscInt component = 0;
                 component < dofPerCell_;
                 ++component)
            {
                indices.push_back(first + component);
            }
        };

    for (PetscInt local = 0;
         local < ownedCellCount_;
         ++local)
    {
        appendCell(firstOwnedCellId_ + local);
    }
    for (PetscInt ghost : ghostIds)
        appendCell(ghost);

    return indices;
}

ISLocalToGlobalMapping
DofMap::createLocalToGlobalMapping() const
{
    const auto blocks =
        localToGlobalBlockIndices();

    ISLocalToGlobalMapping mapping =
        nullptr;

    PetscCallAbort(
        mesh_.communicator(),
        ISLocalToGlobalMappingCreate(
            mesh_.communicator(),
            dofPerCell_,
            static_cast<PetscInt>(
                blocks.size()),
            blocks.data(),
            PETSC_COPY_VALUES,
            &mapping));

    return mapping;
}

} // namespace MPMC
