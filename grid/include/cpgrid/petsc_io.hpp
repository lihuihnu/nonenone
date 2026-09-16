/**
 * @file petsc_io.hpp
 * @brief grid 模块中的 `petsc_io` 源码。
 */
#pragma once

#include <petscmat.h>
#include <petscvec.h>

#include <string>
#include <vector>

namespace MPMC
{

/**
 * @brief 将 PETSc Vec 保存为 MATLAB 兼容 ASCII。
 *
 * 默认变量名与文件名主干一致；可通过 @p variableName 显式覆盖。
 */
void saveVectorMatlab(
    Vec vector,
    const std::string &filename,
    const std::string &variableName = {});

/**
 * @brief 将 PETSc Mat 保存为 MATLAB 兼容 ASCII。
 *
 * 默认变量名与文件名主干一致；可通过 @p variableName 显式覆盖。
 */
void saveMatrixMatlab(
    Mat matrix,
    const std::string &filename,
    const std::string &variableName = {});

/**
 * @brief 将分布式 Vec 收集到指定 root。
 *
 * 非 root rank 返回空 vector。
 */
[[nodiscard]] std::vector<PetscScalar>
gatherVectorToRoot(
    Vec vector,
    PetscMPIInt root = 0);

} // namespace MPMC
