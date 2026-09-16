/**
 * @file dof_layout.hpp
 * @brief CpGrid 单元自由度布局及 PETSc owned/ghost 映射。
 */
#pragma once

#include <cpgrid/dof_map.hpp>
#include <cpgrid/sparsity.hpp>

#include <petscdm.h>
#include <petscis.h>
#include <petscmat.h>
#include <petscsf.h>
#include <petscvec.h>

#include <memory>
#include <vector>

namespace MPMC
{

/**
 * @brief 一个固定每单元 DOF 数的 PETSc 代数布局。
 *
 * DofLayout 使用 PETSc 公共 DMShell 管理：
 * - 全局分布式 Vec；
 * - 本地 `[owned][ghost]` 顺序 Vec；
 * - global-to-local / local-to-global scatter；
 * - block sparse Jacobian。
 *
 * 不访问 PETSc 私有 `dm->ops`、`dm->data`，也不维护自定义 DM 类型。
 */
class DofLayout final
{
  public:
    DofLayout(
        Mesh &mesh,
        PetscInt dofPerCell);

    DofLayout(const DofLayout &) = delete;
    DofLayout &operator=(const DofLayout &) = delete;
    DofLayout(DofLayout &&) = delete;
    DofLayout &operator=(DofLayout &&) = delete;

    ~DofLayout() noexcept;

    [[nodiscard]] PetscInt
    dofPerCell() const noexcept
    {
        return dofMap_.dofPerCell();
    }

    [[nodiscard]] const DofMap &
    dofMap() const noexcept
    {
        return dofMap_;
    }

    [[nodiscard]] DofMap &
    dofMap() noexcept
    {
        return dofMap_;
    }

    [[nodiscard]] DM dm() const noexcept
    {
        return dm_;
    }

    [[nodiscard]] Vec createGlobalVector() const;
    [[nodiscard]] Vec createLocalVector() const;

    [[nodiscard]] Vec createLocalVector(
        Vec global) const;

    /** @brief 从 PETSc DM 临时向量池借用 owned+ghost 局部 Vec。 */
    [[nodiscard]] Vec borrowLocalVector() const;

    /** @brief 借用局部 Vec，并立即完成 global -> local scatter。 */
    [[nodiscard]] Vec borrowLocalVector(Vec global) const;

    /** @brief 将 borrowLocalVector() 获得的 Vec 归还 PETSc DM 向量池。 */
    void restoreLocalVector(Vec &local) const noexcept;

    void globalToLocal(
        Vec global,
        Vec local,
        InsertMode mode = INSERT_VALUES) const;

    void localToGlobal(
        Vec local,
        Vec global,
        InsertMode mode = ADD_VALUES) const;

    [[nodiscard]] Mat createMatrix() const;

  private:
    static PetscErrorCode
    createGlobalVectorCallback_(
        DM dm,
        Vec *vector);

    static PetscErrorCode
    createLocalVectorCallback_(
        DM dm,
        Vec *vector);

    static PetscErrorCode
    createMatrixCallback_(
        DM dm,
        Mat *matrix);

    static PetscErrorCode
    context_(
        DM dm,
        DofLayout **layout);

    PetscErrorCode
    createGlobalVectorImpl_(
        Vec *vector) const;

    PetscErrorCode
    createLocalVectorImpl_(
        Vec *vector) const;

    PetscErrorCode
    createMatrixImpl_(
        Mat *matrix) const;

    void createShell_();
    void createMapping_(const std::vector<PetscInt> &globalIndexValues);
    [[nodiscard]] VecScatter createScatter_(
        bool localToGlobal,
        const std::vector<PetscInt> &globalIndexValues) const;
    void createGlobalToLocalScatter_(
        const std::vector<PetscInt> &globalIndexValues);
    void ensureLocalToGlobalScatter_() const;
    void release_() noexcept;

    Mesh &mesh_;
    DofMap dofMap_;

    // Jacobian 稀疏模式只在真正请求 Matrix 时构造；1/3-DOF 岩石布局
    // 以及 phase-state 布局不再为未使用的 Jacobian 扫描 topology。
    mutable std::unique_ptr<Sparsity> sparsity_;

    DM dm_{nullptr};
    ISLocalToGlobalMapping localToGlobalMapping_{nullptr};
    VecScatter globalToLocalScatter_{nullptr};

    // production CpGrid 只消费 global -> local；反向 scatter 保留旧 API，
    // 第一次 localToGlobal() 时才创建。
    mutable VecScatter localToGlobalScatter_{nullptr};
};

} // namespace MPMC
