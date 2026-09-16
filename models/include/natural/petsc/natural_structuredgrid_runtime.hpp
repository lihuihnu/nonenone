/**
 * @file natural_structuredgrid_runtime.hpp
 * @brief Natural/PETSc 在 StructuredGrid 上的便捷运行时入口。
 */
#pragma once

#include <natural/petsc/structuredgrid_backend.hpp>
#include <natural/petsc/natural_petsc_runtime.hpp>
#include <structuredgrid/structuredgrid.hpp>

namespace MPMC
{

/**
 * @brief Natural + StructuredGrid(DMDA) 的标准运行时入口。
 *
 * StructuredGrid 只提供规则网格、DMDA/ghost、岩石属性和几何；Natural 的
 * EOS/flash/通量/蓄积/井/Jacobian 与 CpGrid 共用 NaturalPetscRuntime。
 */
template <
    class Indices,
    class Grid = StructuredGridCore>
using NaturalStructuredGridRuntime =
    NaturalPetscRuntime<
        Indices,
        NaturalStructuredGridBackend<Indices, Grid>>;

} // namespace MPMC
