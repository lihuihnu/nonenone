/**
 * @file petsc_custom_hooks.hpp
 * @brief 统一导出定制 PETSc 所需的 Natural update hook 及 C/C++ ABI 分发入口。
 *
 * 目标超算上的定制 PETSc 可能以 C++ 语言构建并引用 mangled C++ `updateState`
 * / `updateSol`；本地 PETSc 也可能以 C 语言构建并引用裸 C 符号。每个 case 通过
 * 本宏始终保留历史 C++ ABI，同时额外导出稳定的 C dispatch 名称。若 machine
 * profile 选择 `PETSC_CUSTOM_HOOK_LINKAGE=both`，case/Makefile 会再链接一个很薄的
 * C ABI bridge，把裸 C update 符号转发到同一 dispatch，两个 ABI 因而执行完全相同
 * 的 Natural runtime 逻辑。
 */
#pragma once

#include <natural/petsc/callbacks.hpp>
#include <petscvec.h>

#define MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(RuntimeType)                         \
    extern "C" void mpmcPetscUpdateStateDispatch(Vec x, void *ctx) noexcept        \
    {                                                                                \
        MPMC::naturalUpdateStateHook<RuntimeType>(x, ctx);                           \
    }                                                                                \
    extern "C" void mpmcPetscUpdateSolDispatch(Vec x, Vec y, void *ctx) noexcept   \
    {                                                                                \
        MPMC::naturalUpdateSolHook<RuntimeType>(x, y, ctx);                          \
    }                                                                                \
    void updateState(Vec x, void *ctx) noexcept                                      \
    {                                                                                \
        mpmcPetscUpdateStateDispatch(x, ctx);                                         \
    }                                                                                \
    void updateSol(Vec x, Vec y, void *ctx) noexcept                                 \
    {                                                                                \
        mpmcPetscUpdateSolDispatch(x, y, ctx);                                        \
    }
