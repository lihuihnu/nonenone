/**
 * @file grid_backend_common.hpp
 * @brief Natural 不同网格后端共享的单元/面访问辅助。
 */
#pragma once

#include <petscsys.h>

namespace MPMC
{

/**
 * @brief Natural runtime 使用的统一单元连接信息。
 */
struct NaturalCellConnection final
{
    PetscInt neighborCellId{-1};
    double transmissibility{0.0};
    double gravityTerm{0.0};
};

/**
 * @brief 把后端的 DOF 查询封装成与 runtime 无关的轻量视图。
 */
template <class Backend>
class NaturalDofMapView final
{
public:
    using CellId = typename Backend::CellId;

    NaturalDofMapView(const Backend &backend, PetscInt dof)
        : backend_(backend), dof_(dof)
    {
    }

    [[nodiscard]] PetscInt globalIndex(CellId cell, PetscInt component) const
    {
        return backend_.globalDof(cell, dof_, component);
    }

    [[nodiscard]] PetscInt localIndex(CellId cell, PetscInt component) const
    {
        return backend_.localIndex(cell, dof_, component);
    }

    [[nodiscard]] PetscInt localBlockIndex(CellId cell) const
    {
        return backend_.localBlockIndex(cell, dof_);
    }

    [[nodiscard]] PetscInt ownedCellCount() const
    {
        return backend_.ownedCellCount(dof_);
    }

    [[nodiscard]] PetscInt ghostCellCount() const
    {
        return backend_.ghostCellCount(dof_);
    }

    [[nodiscard]] PetscInt ownedDofCount() const
    {
        return backend_.ownedDofCount(dof_);
    }

    [[nodiscard]] PetscInt globalDofCount() const
    {
        return backend_.globalDofCount(dof_);
    }

private:
    const Backend &backend_;
    PetscInt dof_{0};
};

} // namespace MPMC
