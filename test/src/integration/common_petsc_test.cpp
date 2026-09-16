/**
 * @file common_petsc_test.cpp
 * @brief 集成测试：验证 `common_petsc_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <common/petsc_io.hpp>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <cmath>
#include <stdexcept>

int main(
    int argc,
    char **argv)
{
    PetscInitialize(
        &argc,
        &argv,
        nullptr,
        nullptr);

    int exitCode = 0;

    try
    {
        Vec reference = nullptr;
        Vec loaded = nullptr;
        Mat matrix = nullptr;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecCreateMPI(
                PETSC_COMM_WORLD,
                2,
                PETSC_DECIDE,
                &reference));

        PetscInt begin = 0;
        PetscInt end = 0;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecGetOwnershipRange(
                reference,
                &begin,
                &end));

        PetscScalar *array = nullptr;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecGetArray(
                reference,
                &array));

        for (PetscInt global = begin;
             global < end;
             ++global)
        {
            array[global - begin] =
                static_cast<PetscScalar>(
                    global + 1);
        }

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecRestoreArray(
                reference,
                &array));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MPMC::petsc::saveVectorAsciiMatlab(
                reference,
                "common_vector.m"));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MPMC::petsc::saveVectorBinary(
                reference,
                "common_vector.bin"));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecDuplicate(
                reference,
                &loaded));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecSet(
                loaded,
                0.0));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MPMC::petsc::loadVectorBinary(
                loaded,
                "common_vector.bin"));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecAXPY(
                loaded,
                -1.0,
                reference));

        PetscReal infinityNorm = 0.0;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecNorm(
                loaded,
                NORM_INFINITY,
                &infinityNorm));

        if (infinityNorm > 1.0e-14)
        {
            throw std::runtime_error(
                "PETSc binary vector round-trip failed.");
        }

        PetscInt globalSize = 0;
        PetscInt localSize = 0;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecGetSize(
                reference,
                &globalSize));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecGetLocalSize(
                reference,
                &localSize));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MatCreateAIJ(
                PETSC_COMM_WORLD,
                localSize,
                localSize,
                globalSize,
                globalSize,
                1,
                nullptr,
                0,
                nullptr,
                &matrix));

        for (PetscInt row = begin;
             row < end;
             ++row)
        {
            const PetscScalar value =
                static_cast<PetscScalar>(
                    row + 1);

            PetscCallAbort(
                PETSC_COMM_WORLD,
                MatSetValue(
                    matrix,
                    row,
                    row,
                    value,
                    INSERT_VALUES));
        }

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MatAssemblyBegin(
                matrix,
                MAT_FINAL_ASSEMBLY));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MatAssemblyEnd(
                matrix,
                MAT_FINAL_ASSEMBLY));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MPMC::petsc::saveMatrixAsciiMatlab(
                matrix,
                "common_matrix.m"));

        PetscPrintf(
            PETSC_COMM_WORLD,
            "Common PETSc I/O validation: PASS "
            "(binary round-trip max error = %.3e)\n",
            static_cast<double>(
                infinityNorm));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            MatDestroy(
                &matrix));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecDestroy(
                &loaded));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecDestroy(
                &reference));
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            PETSC_COMM_WORLD,
            "Common PETSc I/O validation failed: %s\n",
            error.what());

        exitCode = 1;
    }

    PetscFinalize();
    return exitCode;
}
