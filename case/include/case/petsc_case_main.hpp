/**
 * @file petsc_case_main.hpp
 * @brief PETSc 储层算例统一 main 生命周期包装。
 */
#pragma once

#include <petscsys.h>

#include <exception>

namespace MPMC::cases
{

/**
 * @brief 统一执行 PETSc 初始化、算例入口、异常报告与 PETSc 收尾。
 *
 * 该包装只收敛各算例重复的进程生命周期，不改变算例 run() 的执行内容。
 * 保持原有行为：仅捕获 std::exception，错误统一打印为 [ERROR][CASE]，
 * 并在正常返回或捕获异常后调用 PetscFinalize()。
 */
template <class RunFunction>
int runPetscCaseMain(int argc, char **argv, RunFunction run)
{
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscInitialize(&argc, &argv, nullptr, nullptr));

    int code = 0;
    try
    {
        code = run();
    }
    catch (const std::exception &error)
    {
        PetscPrintf(PETSC_COMM_WORLD, "[ERROR][CASE] %s\n", error.what());
        code = 1;
    }

    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
    return code;
}

} // namespace MPMC::cases
