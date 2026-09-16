/**
 * @file petsc_custom_hooks_c_bridge.cpp
 * @brief 为以 C 语言构建的定制 PETSc 提供裸 `updateState/updateSol` ABI 转发。
 *
 * 该文件只在 machine profile 选择 `PETSC_CUSTOM_HOOK_LINKAGE=both` 时链接。
 * case 主翻译单元仍保留历史 C++ update 符号；这里的 C 符号调用独立 dispatch，
 * 因而不会递归，也不会复制 Natural 状态更新/步长限制逻辑。
 */
#include <petscvec.h>

extern "C" void mpmcPetscUpdateStateDispatch(Vec, void *) noexcept;
extern "C" void mpmcPetscUpdateSolDispatch(Vec, Vec, void *) noexcept;

extern "C" void updateState(Vec x, void *ctx) noexcept
{
    mpmcPetscUpdateStateDispatch(x, ctx);
}

extern "C" void updateSol(Vec x, Vec y, void *ctx) noexcept
{
    mpmcPetscUpdateSolDispatch(x, y, ctx);
}
