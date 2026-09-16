/**
 * @file structuredgrid_tag.hpp
 * @brief 旧 StructuredGrid<Tag> 兼容 façade 使用的 Natural 编译期标签。
 */
#pragma once

namespace MPMC
{

/**
 * @brief 仅用于兼容旧 `StructuredGrid<Tag>` 源码的 Natural Tag。
 *
 * 新生产路径使用 `StructuredGridCore`，布局由 Natural backend 显式注册；
 * 本类型不再是 Grid 核心依赖。
 */
template <class Indices>
struct NaturalStructuredGridTag final
{
    using ValueType = typename Indices::ValueType;

    static constexpr int numVars_ =
        Indices::numPrimaryVariables;

    static constexpr int numPhaseState_ =
        Indices::numPhaseStateVariables;
};

} // namespace MPMC
