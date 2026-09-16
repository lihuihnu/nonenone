/**
 * @file petsc_io.cpp
 * @brief CpGrid PETSc 数据导入导出的实现。
 */
#include <cpgrid/petsc_io.hpp>
#include <common/petsc_io.hpp>

#include <petscsys.h>
#include <petscsf.h>
#include <petscviewer.h>

#include <stdexcept>

namespace MPMC
{

void saveVectorMatlab(
    Vec vector,
    const std::string &filename,
    const std::string &variableName)
{
    PetscCallAbort(
        vector != nullptr
            ? PetscObjectComm(reinterpret_cast<PetscObject>(vector))
            : PETSC_COMM_SELF,
        MPMC::petsc::saveVectorAsciiMatlab(vector, filename, variableName));
}

void saveMatrixMatlab(
    Mat matrix,
    const std::string &filename,
    const std::string &variableName)
{
    PetscCallAbort(
        matrix != nullptr
            ? PetscObjectComm(reinterpret_cast<PetscObject>(matrix))
            : PETSC_COMM_SELF,
        MPMC::petsc::saveMatrixAsciiMatlab(matrix, filename, variableName));
}

std::vector<PetscScalar>
gatherVectorToRoot(
    Vec vector,
    PetscMPIInt root)
{
    if (vector == nullptr)
    {
        throw std::invalid_argument(
            "gatherVectorToRoot requires a valid Vec.");
    }

    MPI_Comm comm =
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                vector));

    PetscMPIInt rank = 0;
    PetscMPIInt size = 0;

    PetscCallMPIAbort(
        comm,
        MPI_Comm_rank(
            comm,
            &rank));

    PetscCallMPIAbort(
        comm,
        MPI_Comm_size(
            comm,
            &size));

    if (root < 0 ||
        root >= size)
    {
        throw std::out_of_range(
            "Requested gather root is outside communicator.");
    }

    /*
     * VecScatterCreateToZero 固定收集到 rank 0。
     * 对非零 root，使用一般 VecScatter 更合适；当前调试接口只允许 root=0，
     * 避免隐藏额外的自定义 scatter 逻辑。
     */
    if (root != 0)
    {
        throw std::invalid_argument(
            "gatherVectorToRoot currently supports root == 0 only.");
    }

    VecScatter scatter = nullptr;
    Vec sequential = nullptr;

    PetscCallAbort(
        comm,
        VecScatterCreateToZero(
            vector,
            &scatter,
            &sequential));

    PetscCallAbort(
        comm,
        VecScatterBegin(
            scatter,
            vector,
            sequential,
            INSERT_VALUES,
            SCATTER_FORWARD));

    PetscCallAbort(
        comm,
        VecScatterEnd(
            scatter,
            vector,
            sequential,
            INSERT_VALUES,
            SCATTER_FORWARD));

    std::vector<PetscScalar> values;

    if (rank == root)
    {
        PetscInt sizeGlobal = 0;

        PetscCallAbort(
            PETSC_COMM_SELF,
            VecGetSize(
                sequential,
                &sizeGlobal));

        values.resize(
            static_cast<std::size_t>(
                sizeGlobal));

        const PetscScalar *array =
            nullptr;

        PetscCallAbort(
            PETSC_COMM_SELF,
            VecGetArrayRead(
                sequential,
                &array));

        for (PetscInt i = 0;
             i < sizeGlobal;
             ++i)
        {
            values[
                static_cast<std::size_t>(i)] =
                array[i];
        }

        PetscCallAbort(
            PETSC_COMM_SELF,
            VecRestoreArrayRead(
                sequential,
                &array));
    }

    PetscCallAbort(
        comm,
        VecScatterDestroy(
            &scatter));

    if (sequential != nullptr)
    {
        PetscCallAbort(
            PETSC_COMM_SELF,
            VecDestroy(
                &sequential));
    }

    return values;
}

} // namespace MPMC
