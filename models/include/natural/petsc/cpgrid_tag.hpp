/**
 * @file cpgrid_tag.hpp
 * @brief 旧 CpGrid<Tag> 兼容 façade 使用的 Natural 编译期标签。
 */
#pragma once

#include <petscsys.h>

namespace MPMC
{

/**
 * @brief 仅用于兼容旧 `CpGrid<Tag>` 源码的主变量 DOF Tag。
 *
 * 新生产路径使用 `CpGridCore`，主变量布局由 Natural backend 显式注册。
 */
template <class Indices>
struct NaturalCpGridTag final
{
    static constexpr PetscInt numVars_ =
        static_cast<PetscInt>(
            Indices::numPrimaryVariables);
};

} // namespace MPMC
