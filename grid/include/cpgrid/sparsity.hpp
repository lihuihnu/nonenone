/**
 * @file sparsity.hpp
 * @brief 由网格邻接关系构造 PETSc Jacobian 稀疏结构。
 */
#pragma once

#include <cpgrid/dof_map.hpp>

#include <petscmat.h>

#include <vector>

namespace MPMC
{

/**
 * @brief CpGrid 单元块 Jacobian 的稀疏结构预计算。
 *
 * 一个 Polyhedron 对应一个 block row；单元自身以及所有内部面邻居都预留一个
 * dofPerCell × dofPerCell 的稠密块。
 */
class Sparsity final
{
  public:
    explicit Sparsity(
        const DofMap &dofMap);

    [[nodiscard]] const std::vector<PetscInt> &
    diagonalBlockNnz() const noexcept
    {
        return diagonalBlockNnz_;
    }

    [[nodiscard]] const std::vector<PetscInt> &
    offDiagonalBlockNnz() const noexcept
    {
        return offDiagonalBlockNnz_;
    }

    /**
     * @brief 创建已经完成块预分配的 PETSc Matrix。
     *
     * 默认类型为 MATBAIJ；用户仍可通过 `-mat_type` 覆盖。
     * 本接口用于 DMShell 的 PETSc 回调，因此使用 PetscErrorCode 向上传递错误，
     * 不在 PETSc callback 内部调用 PetscCallAbort()。
     */
    PetscErrorCode createMatrix(
        Mat *matrix) const;

  private:
    void build_();

    const DofMap &dofMap_;
    std::vector<PetscInt> diagonalBlockNnz_;
    std::vector<PetscInt> offDiagonalBlockNnz_;
};

} // namespace MPMC
