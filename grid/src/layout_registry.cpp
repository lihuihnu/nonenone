/**
 * @file layout_registry.cpp
 * @brief 不同每单元自由度布局的 PETSc 对象缓存与生命周期管理的实现。
 */
#include <cpgrid/layout_registry.hpp>

#include <stdexcept>

namespace MPMC
{

DofLayout &LayoutRegistry::registerLayout(
    PetscInt dofPerCell)
{
    if (dofPerCell <= 0)
    {
        throw std::invalid_argument(
            "LayoutRegistry requires dofPerCell > 0.");
    }

    const auto found =
        layouts_.find(dofPerCell);

    if (found != layouts_.end())
    {
        return *found->second;
    }

    auto layout =
        std::make_unique<DofLayout>(
            mesh_,
            dofPerCell);

    DofLayout &reference = *layout;

    layouts_.emplace(
        dofPerCell,
        std::move(layout));

    return reference;
}

bool LayoutRegistry::contains(
    PetscInt dofPerCell) const noexcept
{
    return layouts_.find(dofPerCell) !=
           layouts_.end();
}

DofLayout &LayoutRegistry::layout(
    PetscInt dofPerCell)
{
    const auto found =
        layouts_.find(dofPerCell);

    if (found == layouts_.end())
    {
        throw std::out_of_range(
            "Requested CpGrid DOF layout has not been registered.");
    }

    return *found->second;
}

const DofLayout &LayoutRegistry::layout(
    PetscInt dofPerCell) const
{
    const auto found =
        layouts_.find(dofPerCell);

    if (found == layouts_.end())
    {
        throw std::out_of_range(
            "Requested CpGrid DOF layout has not been registered.");
    }

    return *found->second;
}

} // namespace MPMC
