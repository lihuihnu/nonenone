/**
 * @file natural_cpgrid_runtime.hpp
 * @brief Natural/PETSc 在 CpGrid 上的便捷运行时别名与入口。
 */
#pragma once

#include <cpgrid/cpgrid.hpp>
#include <natural/petsc/cpgrid_backend.hpp>
#include <natural/petsc/natural_petsc_runtime.hpp>

namespace MPMC
{

/**
 * @brief Natural + CpGrid 的标准入口。
 *
 * 真正的 residual/Jacobian/updateState/updateSol 实现在 NaturalPetscRuntime；
 * 此别名只选择 CpGrid backend，为用户隐藏内部 backend 与数值常量类型。
 */
template <
    class Indices,
    class Grid = CpGridCore>
using NaturalCpGridRuntime =
    NaturalPetscRuntime<
        Indices,
        NaturalCpGridBackend<Indices, Grid>>;

} // namespace MPMC
