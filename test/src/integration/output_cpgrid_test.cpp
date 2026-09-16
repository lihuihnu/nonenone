/**
 * @file output_cpgrid_test.cpp
 * @brief 集成测试：验证 `output_cpgrid_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <output/petsc/output_petsc.hpp>

#include <cpgrid/cpgrid.hpp>
#include <cpgrid/distributed_mesh_loader.hpp>
#include <indices/indices.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Config =
    MPMC::CompositionalModelConfig<
        6,
        true,
        true,
        false,
        false,
        false>;

using Indices =
    MPMC::ScalarIndices<Config>;

using Grid = MPMC::CpGridCore;

double fieldValue(
    double base,
    PetscInt currentCell,
    PetscInt component)
{
    return base +
           10.0 *
               static_cast<double>(
                   currentCell) +
           static_cast<double>(
               component);
}

void fillField(
    Grid &grid,
    Vec vector,
    PetscInt dofPerCell,
    double base)
{
    auto &map =
        grid.dofMap(
            dofPerCell);

    PetscScalar *array = nullptr;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetArray(
            vector,
            &array));

    for (const auto &cell :
         grid.localCells())
    {
        const PetscInt localBlock =
            map.localBlockIndex(
                cell);

        for (PetscInt component = 0;
             component < dofPerCell;
             ++component)
        {
            array[
                localBlock *
                    dofPerCell +
                component] =
                static_cast<PetscScalar>(
                    fieldValue(
                        base,
                        cell.id(),
                        component));
        }
    }

    PetscCallAbort(
        grid.mesh().communicator(),
        VecRestoreArray(
            vector,
            &array));
}

double componentTotalReference(
    PetscInt cellCount,
    PetscInt component,
    double base)
{
    const double n =
        static_cast<double>(
            cellCount);

    return
        n *
            (base +
             static_cast<double>(
                 component)) +
        10.0 *
            n *
            (n - 1.0) /
            2.0;
}

void run(
    const std::string &meshDirectory)
{
    auto distributed = MPMC::loadDistributedMeshFromRoot(
        meshDirectory,
        false,
        {},
        PETSC_COMM_WORLD,
        0);
    MPMC::Mesh &mesh = *distributed.mesh;

    // 输出采样只依赖全局轻量编号目录，不要求 root 长期保留远端 Polyhedron。
    mesh.compactToLocalSnapshot(false);

    Grid grid(mesh);

    grid.registerLayout(
        static_cast<PetscInt>(
            Indices::numPrimaryVariables));

    grid.registerLayout(
        static_cast<PetscInt>(
            Indices::numPhaseStateVariables));

    grid.registerLayout(
        static_cast<PetscInt>(
            Indices::numComponents));

    Vec solution =
        grid.createGlobalVector(
            static_cast<PetscInt>(
                Indices::numPrimaryVariables));

    Vec phaseState =
        grid.createGlobalVector(
            static_cast<PetscInt>(
                Indices::numPhaseStateVariables));

    Vec componentMass =
        grid.createGlobalVector(
            static_cast<PetscInt>(
                Indices::numComponents));

    fillField(
        grid,
        solution,
        static_cast<PetscInt>(
            Indices::numPrimaryVariables),
        100.0);

    fillField(
        grid,
        phaseState,
        static_cast<PetscInt>(
            Indices::numPhaseStateVariables),
        200.0);

    fillField(
        grid,
        componentMass,
        static_cast<PetscInt>(
            Indices::numComponents),
        300.0);

    MPMC::CpGridOutputAccess<Grid>
        access(grid);

    const PetscInt cellCount =
        static_cast<PetscInt>(
            mesh.cellCount());

    const std::vector<
        MPMC::OutputCellSelection>
        selections{
            {0, MPMC::OutputCellIdSpace::CurrentId},
            {cellCount - 1, MPMC::OutputCellIdSpace::CurrentId},
            {0, MPMC::OutputCellIdSpace::InputIndex}};

    const std::vector<PetscInt>
        variables{
            0,
            static_cast<PetscInt>(
                Indices::numPrimaryVariables - 1)};

    const auto sampled =
        access.sample(
            solution,
            static_cast<PetscInt>(
                Indices::numPrimaryVariables),
            selections,
            variables);

    const PetscInt inputZeroCurrent =
        mesh.currentIdFromInputIndex(0);

    const std::array<PetscInt, 3>
        expectedCells{
            0,
            cellCount - 1,
            inputZeroCurrent};

    double sampleError = 0.0;

    for (std::size_t cell = 0;
         cell < expectedCells.size();
         ++cell)
    {
        for (std::size_t variable = 0;
             variable < variables.size();
             ++variable)
        {
            const double expected =
                fieldValue(
                    100.0,
                    expectedCells[cell],
                    variables[variable]);

            const double actual =
                sampled[
                    cell *
                        variables.size() +
                    variable];

            sampleError =
                std::max(
                    sampleError,
                    std::abs(
                        actual -
                        expected));
        }
    }

    // The solver vector is current-id ordered internally.  External copies must
    // be input-row ordered regardless of the METIS partition.
    Vec inputOrderedSolution =
        grid.createInputOrderedCopy(
            static_cast<PetscInt>(
                Indices::numPrimaryVariables),
            solution);

    Vec gatheredInputOrdered = nullptr;
    VecScatter gatherScatter = nullptr;
    PetscCallAbort(
        mesh.communicator(),
        VecScatterCreateToZero(
            inputOrderedSolution,
            &gatherScatter,
            &gatheredInputOrdered));
    PetscCallAbort(
        mesh.communicator(),
        VecScatterBegin(
            gatherScatter,
            inputOrderedSolution,
            gatheredInputOrdered,
            INSERT_VALUES,
            SCATTER_FORWARD));
    PetscCallAbort(
        mesh.communicator(),
        VecScatterEnd(
            gatherScatter,
            inputOrderedSolution,
            gatheredInputOrdered,
            INSERT_VALUES,
            SCATTER_FORWARD));

    double inputOrderError = 0.0;
    if (mesh.rank() == 0)
    {
        const PetscScalar *ordered = nullptr;
        PetscCallAbort(
            PETSC_COMM_SELF,
            VecGetArrayRead(
                gatheredInputOrdered,
                &ordered));

        for (PetscInt inputCell = 0; inputCell < cellCount; ++inputCell)
        {
            const PetscInt currentCell =
                mesh.currentIdFromInputIndex(
                    static_cast<std::size_t>(inputCell));

            for (PetscInt component = 0;
                 component < Indices::numPrimaryVariables;
                 ++component)
            {
                const PetscInt global =
                    inputCell *
                        static_cast<PetscInt>(Indices::numPrimaryVariables) +
                    component;

                inputOrderError = std::max(
                    inputOrderError,
                    std::abs(
                        PetscRealPart(ordered[global]) -
                        fieldValue(100.0, currentCell, component)));
            }
        }

        PetscCallAbort(
            PETSC_COMM_SELF,
            VecRestoreArrayRead(
                gatheredInputOrdered,
                &ordered));
    }

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Bcast(
            &inputOrderError,
            1,
            MPI_DOUBLE,
            0,
            mesh.communicator()));

    PetscCallAbort(
        mesh.communicator(),
        VecScatterDestroy(&gatherScatter));
    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(&gatheredInputOrdered));
    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(&inputOrderedSolution));

    const auto totals =
        access.template componentTotals<
            static_cast<std::size_t>(
                Indices::numComponents)>(
                    componentMass);

    double totalError = 0.0;

    for (PetscInt component = 0;
         component < Indices::numComponents;
         ++component)
    {
        const double expected =
            componentTotalReference(
                cellCount,
                component,
                300.0);

        totalError =
            std::max(
                totalError,
                std::abs(
                    totals.component[
                        static_cast<std::size_t>(
                            component)] -
                    expected));
    }

    {
        MPMC::CpGridSaver<
            Indices,
            Grid>
            saver(
                "output_cpgrid_test_results",
                grid);

        saver.saveSolution(
            solution,
            {
                {0, MPMC::OutputCellIdSpace::CurrentId},
                {cellCount - 1, MPMC::OutputCellIdSpace::CurrentId}},
            {
                0,
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables - 1)});

        saver.savePhaseState(
            phaseState,
            {
                {0, MPMC::OutputCellIdSpace::CurrentId}},
            {
                static_cast<PetscInt>(
                    Indices::PhaseState::flag),
                static_cast<PetscInt>(
                    Indices::PhaseState::liquidFraction)});

        saver.saveMass(
            MPMC::MassSeries::OilPhaseComponents,
            componentMass,
            {
                {0, MPMC::OutputCellIdSpace::CurrentId}},
            {0, 1});

        saver.saveMassTotal(
            MPMC::MassSeries::OilPhaseComponents,
            componentMass,
            2.5);
    }

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Barrier(
            mesh.communicator()));

    if (mesh.rank() == 0)
    {
        if (!std::filesystem::exists(
                "output_cpgrid_test_results/solution_samples.csv") ||
            !std::filesystem::exists(
                "output_cpgrid_test_results/oil_phase_component_mass_totals.csv"))
        {
            throw std::runtime_error(
                "CpGrid output writer did not create expected rank-0 files.");
        }
    }

    if (sampleError > 1.0e-12 ||
        inputOrderError > 1.0e-12 ||
        totalError > 1.0e-8)
    {
        throw std::runtime_error(
            "CpGrid output MPI value validation failed.");
    }

    PetscPrintf(
        mesh.communicator(),
        "Output/CpGrid MPI validation\n"
        "  ranks             = %d\n"
        "  cells             = %" PetscInt_FMT "\n"
        "  sample max error  = %.3e\n"
        "  input-order error = %.3e\n"
        "  total max error   = %.3e\n"
        "  result            = ALL PASS\n",
        mesh.processCount(),
        cellCount,
        sampleError,
        inputOrderError,
        totalError);

    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(
            &componentMass));

    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(
            &phaseState));

    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(
            &solution));

    if (mesh.rank() == 0)
    {
        std::filesystem::remove_all(
            "output_cpgrid_test_results");
    }
}

} // namespace

int main(
    int argc,
    char **argv)
{
    PetscInitialize(
        &argc,
        &argv,
        nullptr,
        nullptr);

    int status = 0;

    try
    {
        char meshDirectory[
            PETSC_MAX_PATH_LEN] = {};

        PetscBool found =
            PETSC_FALSE;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetString(
                nullptr,
                nullptr,
                "-mesh_dir",
                meshDirectory,
                sizeof(meshDirectory),
                &found));

        if (!found)
        {
            throw std::invalid_argument(
                "Please provide -mesh_dir <MRST-export-directory>.");
        }

        {
            run(
                std::string(
                    meshDirectory));
        }
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            PETSC_COMM_WORLD,
            "Output/CpGrid test failed: %s\n",
            error.what());

        status = 1;
    }

    PetscFinalize();
    return status;
}
