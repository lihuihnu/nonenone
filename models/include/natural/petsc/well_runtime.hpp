/**
 * @file well_runtime.hpp
 * @brief Natural 运行时中的井状态更新、控制检查与源项耦合。
 */
#pragma once

#include <natural/physics/well_source.hpp>
#include <well/specification.hpp>

#include <petscsys.h>

namespace MPMC
{

/**
 * @brief Natural PETSc runtime 使用的 canonical well 类型。
 *
 * 井数据模型与网格无关：StructuredGrid 使用 Cartesian cell id，CpGrid 使用
 * current cell id；二者都以 PetscInt 存储，由对应 grid backend 解释。
 */
template <class Indices>
using NaturalPerforation = WellPerforation<PetscInt>;

template <class Indices>
using NaturalWell = WellSpecification<Indices, PetscInt>;


} // namespace MPMC
