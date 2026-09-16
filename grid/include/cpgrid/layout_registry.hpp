/**
 * @file layout_registry.hpp
 * @brief 不同每单元自由度布局的 PETSc 对象缓存与生命周期管理。
 */
#pragma once

#include <cpgrid/dof_layout.hpp>

#include <petscsys.h>

#include <cstddef>
#include <memory>
#include <unordered_map>

namespace MPMC
{

/**
 * @brief 管理同一 Mesh 上不同每单元 DOF 数的 DofLayout。
 *
 * DMShell 创建以及 scatter 构造包含 MPI 集体操作，因此所有 rank 必须以一致
 * 的顺序注册相同 DOF 布局。
 */
class LayoutRegistry final
{
  public:
    explicit LayoutRegistry(
        Mesh &mesh) noexcept
        : mesh_(mesh)
    {
    }

    LayoutRegistry(
        const LayoutRegistry &) = delete;
    LayoutRegistry &operator=(
        const LayoutRegistry &) = delete;
    LayoutRegistry(
        LayoutRegistry &&) = delete;
    LayoutRegistry &operator=(
        LayoutRegistry &&) = delete;

    ~LayoutRegistry() = default;

    DofLayout &registerLayout(
        PetscInt dofPerCell);

    [[nodiscard]] bool contains(
        PetscInt dofPerCell) const noexcept;

    [[nodiscard]] DofLayout &layout(
        PetscInt dofPerCell);

    [[nodiscard]] const DofLayout &layout(
        PetscInt dofPerCell) const;

    [[nodiscard]] std::size_t size() const noexcept
    {
        return layouts_.size();
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return layouts_.empty();
    }

  private:
    Mesh &mesh_;
    std::unordered_map<
        PetscInt,
        std::unique_ptr<DofLayout>>
        layouts_;
};

} // namespace MPMC
