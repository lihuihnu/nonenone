/**
 * @file petsc_io.hpp
 * @brief common 模块中的 `petsc_io` 源码。
 */
#pragma once

#include <petscmat.h>
#include <petscvec.h>

#include <string>
#include <vector>

namespace MPMC::petsc
{

/**
 * @brief 以 MATLAB ASCII 格式保存 PETSc Vec。
 *
 * 使用 Vec 自身的 communicator，不假定 PETSC_COMM_WORLD。
 *
 * MATLAB 变量名默认由文件名自动生成：去掉目录和扩展名，并规范为合法
 * MATLAB 标识符。例如 `solution_step_60.m` 内保存变量
 * `solution_step_60`。如果传入 @p variableName，则使用显式变量名。
 */
PetscErrorCode saveVectorAsciiMatlab(
    Vec vector,
    const std::string &filename,
    const std::string &variableName = {});

/**
 * @brief 以 MATLAB ASCII 格式保存 PETSc Mat。
 *
 * 使用 Mat 自身的 communicator，不假定 PETSC_COMM_WORLD。
 *
 * 默认变量名同样与文件名主干一致；也可通过 @p variableName 显式指定。
 */
PetscErrorCode saveMatrixAsciiMatlab(
    Mat matrix,
    const std::string &filename,
    const std::string &variableName = {});


/**
 * @brief 将按单元分块的 PETSc Vec 保存为常规 CSV 表。
 *
 * The input Vec is interpreted as contiguous cell blocks:
 * `cell0[dof0..dofN-1], cell1[dof0..dofN-1], ...`.  The first
 * CSV column is the zero-based cell index, followed by one column per DOF.
 * This routine does not reorder cells; callers that solve in a partition-dependent
 * ordering must first create an input-ordered copy.
 *
 * @param vector PETSc vector to save.
 * @param dofPerCell Number of consecutive values stored for each cell.
 * @param filename Output CSV path.
 * @param columnNames Optional names for the DOF columns. If empty, `value_0`,
 *        `value_1`, ... are generated. When provided, its size must equal
 *        @p dofPerCell.
 */
PetscErrorCode saveCellVectorCsv(
    Vec vector,
    PetscInt dofPerCell,
    const std::string &filename,
    const std::vector<std::string> &columnNames = {});

/**
 * @brief 以 PETSc binary 格式保存 Vec。
 */
PetscErrorCode saveVectorBinary(
    Vec vector,
    const std::string &filename);

/**
 * @brief 从 PETSc binary 文件加载到已经创建好的 Vec。
 */
PetscErrorCode loadVectorBinary(
    Vec vector,
    const std::string &filename);

} // namespace MPMC::petsc
