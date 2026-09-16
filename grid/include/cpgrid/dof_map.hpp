/**
 * @file dof_map.hpp
 * @brief 分区网格单元与 PETSc 全局/局部自由度编号映射。
 */
#pragma once

#include <cpgrid/polyhedron.hpp>
#include <cpgrid/mesh.hpp>

#include <petscis.h>
#include <petscsys.h>

#include <vector>

namespace MPMC
{

/**
 * @brief 固定“每单元 DOF 数”的并行自由度映射。
 *
 * 全局排列固定为：
 *
 * ```text
 * rank -> cell -> component
 * ```
 *
 * 本地快照排列固定为：
 *
 * ```text
 * [owned cells][ghost cells]
 * ```
 *
 * 每个 cell 内始终保留 dofPerCell 个连续 DOF，因此与 BAIJ block matrix 和
 * `MatSetValuesBlockedLocal()` 的块编号保持一致。
 */
class DofMap final
{
  public:
    DofMap(
        Mesh &mesh,
        PetscInt dofPerCell);

    DofMap(const DofMap &) = delete;
    DofMap &operator=(const DofMap &) = delete;
    DofMap(DofMap &&) = delete;
    DofMap &operator=(DofMap &&) = delete;

    ~DofMap() = default;

    [[nodiscard]] Mesh &mesh() noexcept
    {
        return mesh_;
    }

    [[nodiscard]] const Mesh &mesh() const noexcept
    {
        return mesh_;
    }

    [[nodiscard]] PetscInt
    dofPerCell() const noexcept
    {
        return dofPerCell_;
    }

    [[nodiscard]] PetscInt
    globalDofCount() const noexcept
    {
        return globalDofCount_;
    }

    [[nodiscard]] PetscInt
    ownedDofCount() const noexcept
    {
        return ownedCellCount_ *
               dofPerCell_;
    }

    [[nodiscard]] PetscInt
    ghostDofCount() const noexcept
    {
        return ghostCellCount() *
               dofPerCell_;
    }

    [[nodiscard]] PetscInt
    localSnapshotDofCount() const noexcept
    {
        return ownedDofCount() +
               ghostDofCount();
    }

    [[nodiscard]] PetscInt
    ownedCellCount() const noexcept
    {
        return ownedCellCount_;
    }

    [[nodiscard]] PetscInt
    ghostCellCount() const noexcept
    {
        return static_cast<PetscInt>(
            mesh_.ghostCellIds().size());
    }

    [[nodiscard]] PetscInt
    firstOwnedCellId() const noexcept
    {
        return firstOwnedCellId_;
    }

    [[nodiscard]] PetscInt
    firstOwnedGlobalDof() const noexcept
    {
        return firstOwnedCellId_ *
               dofPerCell_;
    }

    [[nodiscard]] PetscInt
    endOwnedGlobalDof() const noexcept
    {
        return firstOwnedGlobalDof() +
               ownedDofCount();
    }

    /**
     * @brief 单元某分量的全局标量 DOF。
     */
    [[nodiscard]] PetscInt globalIndex(
        const Polyhedron &cell,
        PetscInt component) const;

    /** @brief 仅凭 current cell id 计算全局标量 DOF，不要求 cell 已物化。 */
    [[nodiscard]] PetscInt globalIndex(
        PetscInt currentCellId,
        PetscInt component) const;

    /**
     * @brief 单元某分量在 `[owned][ghost]` 本地快照中的标量下标。
     */
    [[nodiscard]] PetscInt localIndex(
        const Polyhedron &cell,
        PetscInt component) const;

    /** @brief current id 在 `[owned][ghost]` 快照中的标量下标。 */
    [[nodiscard]] PetscInt localIndex(
        PetscInt currentCellId,
        PetscInt component) const;

    /**
     * @brief 单元的全局 block id，即 current cell id。
     */
    [[nodiscard]] PetscInt globalBlockIndex(
        const Polyhedron &cell) const noexcept
    {
        return cell.id();
    }

    /**
     * @brief 单元在本地 `[owned][ghost]` cell block 排列中的 block id。
     */
    [[nodiscard]] PetscInt localBlockIndex(
        const Polyhedron &cell) const;

    /** @brief current id 在本地 `[owned][ghost]` cell block 中的下标。 */
    [[nodiscard]] PetscInt localBlockIndex(
        PetscInt currentCellId) const;

    [[nodiscard]] bool isOwned(
        const Polyhedron &cell) const noexcept
    {
        return isOwnedCellId(cell.id());
    }

    /** @brief 判断 current cell id 是否落在当前 rank 连续 owned 区间。 */
    [[nodiscard]] bool isOwnedCellId(PetscInt currentCellId) const noexcept
    {
        return currentCellId >= firstOwnedCellId_ &&
               currentCellId < firstOwnedCellId_ + ownedCellCount_;
    }

    [[nodiscard]] bool isOwnedGlobalDof(
        PetscInt globalDof) const noexcept
    {
        return globalDof >=
                   firstOwnedGlobalDof() &&
               globalDof <
                   endOwnedGlobalDof();
    }

    /**
     * @brief 当前 rank 所需的远端 current cell ids。
     */
    [[nodiscard]] const std::vector<PetscInt> &
    ghostCellIds() const noexcept
    {
        return mesh_.ghostCellIds();
    }

    /**
     * @brief 本地 cell block -> 全局 cell block 的映射。
     *
     * 结果为 `[owned current cell ids][ghost current cell ids]`。
     */
    [[nodiscard]] std::vector<PetscInt>
    localToGlobalBlockIndices() const;

    /**
     * @brief 本地标量 DOF -> 全局标量 DOF 的映射。
     */
    [[nodiscard]] std::vector<PetscInt>
    localToGlobalIndices() const;

    /**
     * @brief 创建 block-size=dofPerCell 的 PETSc local-to-global mapping。
     *
     * 返回对象由调用方使用 `ISLocalToGlobalMappingDestroy()` 释放。
     */
    [[nodiscard]] ISLocalToGlobalMapping
    createLocalToGlobalMapping() const;

  private:
    void buildOwnership_();
    void buildGhostMap_();
    void validateComponent_(
        PetscInt component) const;

    Mesh &mesh_;
    PetscInt dofPerCell_{0};

    PetscInt globalDofCount_{0};
    PetscInt ownedCellCount_{0};
    PetscInt firstOwnedCellId_{0};
};

} // namespace MPMC
