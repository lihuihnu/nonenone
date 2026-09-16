#include <petscvec.h>

/**
 * @file petsc_custom_hooks.cpp
 * @brief 测试支撑代码：为 PETSc 测试提供 `petsc_custom_hooks` 相关钩子。
 *
 * 当前超算上的 PETSc 3.22.2 在 SNES 内部直接引用：
 *
 *   void updateState(Vec, void *);
 *   void updateSol(Vec, Vec, void *);
 *
 * 因此，即使某个测试完全不使用 Natural/SNES，链接 libpetsc.so 时也必须能找到
 * 这两个符号。这里提供“测试专用”的弱空实现，仅用于让普通 PETSc/Grid/Output
 * 测试能够链接。
 *
 * Natural + CpGrid 的集成测试会在自己的翻译单元中提供同名强实现，并将调用真正
 * 转发给 Natural runtime。弱符号会自动被强符号覆盖，不改变 Natural 的实际行为。
 *
 * 重要：这些函数属于当前定制 PETSc 的测试环境 ABI，不属于 MPMC_SCW 公共接口，
 * 因而只放在 test/ 中，绝不能移入 common/models 等生产模块。
 */

#if defined(__GNUC__) || defined(__clang__)
#  define MPMC_TEST_WEAK __attribute__((weak))
#else
#  define MPMC_TEST_WEAK
#endif

MPMC_TEST_WEAK void updateState(Vec, void *)
{
    // 普通测试不使用定制 SNES 的 Natural 状态更新；故意为空。
}

MPMC_TEST_WEAK void updateSol(Vec, Vec, void *)
{
    // 普通测试不使用定制 SNES 的 Natural Newton-step limiter；故意为空。
}

#undef MPMC_TEST_WEAK
