/**
 * @file natural_layout_registration.hpp
 * @brief Natural 模型向网格 backend 显式注册代数布局的轻量策略。
 */
#pragma once

#include <type_traits>
#include <utility>

namespace MPMC
{
namespace detail
{

template <class Backend, class = void>
struct HasAdAuxiliaryLayoutRegistration : std::false_type
{
};

template <class Backend>
struct HasAdAuxiliaryLayoutRegistration<
    Backend,
    std::void_t<decltype(
        std::declval<Backend &>().registerAdAuxiliaryLayout(0))>>
    : std::true_type
{
};

template <class Backend>
void registerAdAuxiliaryLayoutIfSupported(Backend &backend, int dof)
{
    if constexpr (HasAdAuxiliaryLayoutRegistration<Backend>::value)
        backend.registerAdAuxiliaryLayout(dof);
}

} // namespace detail

/**
 * @brief 将 Natural 自身需要的每单元 DOF 布局注册到 backend。
 *
 * GridCore 不再读取 Natural Tag。主未知量和 phase-state 由 Natural 层显式
 * 注册；StructuredGrid 历史 AD 辅助布局通过可选 backend hook 保留。没有该
 * hook 的第三方/旧 backend 仍保持原接口兼容。
 *
 * Backend 必须提供 `registerLayout(dof)`；可选提供
 * `registerAdAuxiliaryLayout(dof)`。
 */
template <class Indices, class Backend>
void registerNaturalGridLayouts(Backend &backend)
{
    backend.registerLayout(Indices::numPrimaryVariables);

    if constexpr (Indices::useAutomaticDifferentiation)
    {
        detail::registerAdAuxiliaryLayoutIfSupported(
            backend,
            1 + Indices::numPrimaryVariables);
    }

    backend.registerLayout(Indices::numPhaseStateVariables);
}

} // namespace MPMC
