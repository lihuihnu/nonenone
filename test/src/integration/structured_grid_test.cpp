/**
 * @file structured_grid_test.cpp
 * @brief 集成测试：验证 `structured_grid_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <structuredgrid/structuredgrid.hpp>
#include <structuredgrid/grid_report.hpp>
#include <test/test_utils.hpp>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <exception>
#include <limits>
#include <utility>
#include <vector>

namespace
{

struct TestTag
{
    using ValueType = double;
    static constexpr int numVars_ = 4;
    static constexpr int numPhaseState_ = 3;
};

using Grid = MPMC::StructuredGrid<TestTag>;

void runTest()
{
    bool rejectedNonFiniteActivity = false;
    try
    {
        Grid invalidGrid(1, 1, 1, MPMC::CellSize{1.0, 1.0, 1.0});
        invalidGrid.setActivityInput(
            {std::numeric_limits<double>::quiet_NaN()});
        invalidGrid.setup();
    }
    catch (const std::invalid_argument &)
    {
        rejectedNonFiniteActivity = true;
    }
    MPMC::Test::require(rejectedNonFiniteActivity,
                        "StructuredGrid accepted a non-finite activity mask");

    const std::vector<double> dx{10.0, 12.0, 8.0, 15.0};
    const double dy = 20.0;
    const std::vector<double> dz{2.0, 3.0};

    Grid grid(4, 3, 2, dx, dy, dz);

    // 24 个单元中关闭 1 个，用于验证 active-cell 处理。
    std::vector<double> activity(24, 1.0);
    activity[grid.ijkToGlobal(1, 1, 0)] = 0.0;
    grid.setActivityInput(std::move(activity));
    grid.setup();

    MPMC::Test::require(grid.isSetup(), "StructuredGrid setup() did not finish");
    const auto dimensions = grid.dimensions();
    MPMC::Test::require(dimensions == std::array<int, 3>{4, 3, 2},
                        "StructuredGrid dimensions changed");
    const auto physicalExtent = grid.physicalExtent();
    MPMC::Test::requireNear(physicalExtent[0], 45.0, 1.0e-12, "Lx mismatch");
    MPMC::Test::requireNear(physicalExtent[1], 60.0, 1.0e-12, "Ly mismatch");
    MPMC::Test::requireNear(physicalExtent[2], 5.0, 1.0e-12, "Lz mismatch");
    MPMC::Test::require(!grid.onProcess(-1) && !grid.onProcess(24),
                        "StructuredGrid onProcess accepted an invalid global index");

    const int global = grid.ijkToGlobal(2, 1, 0);
    const auto ijk = grid.globalToIJK(global);
    MPMC::Test::require(ijk == std::array<int, 3>{2, 1, 0},
                        "StructuredGrid index round-trip failed");

    const auto cellSize = grid.cellSize({1, 2, 1});
    MPMC::Test::requireNear(cellSize[0], 12.0, 1.0e-12, "dx mismatch");
    MPMC::Test::requireNear(cellSize[1], 20.0, 1.0e-12, "dy mismatch");
    MPMC::Test::requireNear(cellSize[2], 3.0, 1.0e-12, "dz mismatch");

    PetscInt localActive = static_cast<PetscInt>(grid.cellIndices().size());
    PetscInt globalActive = 0;
    PetscCallMPIAbort(PETSC_COMM_WORLD,
                      MPI_Allreduce(&localActive,
                                    &globalActive,
                                    1,
                                    MPIU_INT,
                                    MPI_SUM,
                                    PETSC_COMM_WORLD));
    MPMC::Test::require(globalActive == 23,
                        "StructuredGrid active-cell count mismatch");

    bool rejectedUnmappedAccess = false;
    try
    {
        (void)grid.cellVolume(global);
    }
    catch (const std::logic_error &)
    {
        rejectedUnmappedAccess = true;
    }
    MPMC::Test::require(rejectedUnmappedAccess,
                        "StructuredGrid allowed geometry access outside an array mapping");

    grid.initializeRockProperties();
    MPMC::Test::require(grid.hasRockProperties(),
                        "StructuredGrid rock-property initialization failed");
    PetscCallAbort(
        grid.communicator(), VecSet(grid.permeabilityVectorHandle(), 2.0));
    PetscCallAbort(
        grid.communicator(), VecSet(grid.porosityVectorHandle(), 0.25));

    {
        auto assemblyAccess = grid.assemblyAccess();
        if (grid.onProcess(global))
        {
            const auto faces = grid.faces(global);
            MPMC::Test::require(faces.size() == 4,
                                "StructuredGrid active face count mismatch");

            const auto xPlus = std::find_if(
                faces.begin(), faces.end(), [](const MPMC::StructuredFace &face) {
                    return face.direction == MPMC::StructuredFaceOrder::XPlus;
                });
            const auto zPlus = std::find_if(
                faces.begin(), faces.end(), [](const MPMC::StructuredFace &face) {
                    return face.direction == MPMC::StructuredFaceOrder::ZPlus;
                });
            MPMC::Test::require(xPlus != faces.end() && zPlus != faces.end(),
                                "StructuredGrid directional faces are incomplete");

            const double expectedTransmissibility =
                40.0 / (4.0 / 2.0 + 7.5 / 2.0);
            MPMC::Test::requireNear(
                grid.transmissibility(global, *xPlus),
                expectedTransmissibility,
                1.0e-12,
                "StructuredGrid directional transmissibility mismatch");
            MPMC::Test::requireNear(
                grid.transmissibility(global, xPlus->neighbor()),
                expectedTransmissibility,
                1.0e-12,
                "StructuredGrid compatibility transmissibility mismatch");
            MPMC::Test::requireNear(
                grid.gravityPotentialDifference(global, *zPlus),
                -2.5 * 9.80665,
                1.0e-12,
                "StructuredGrid directional gravity mismatch");
            MPMC::Test::requireNear(
                grid.cellVolume(global), 320.0, 1.0e-12,
                "StructuredGrid cell volume mismatch");
            MPMC::Test::requireNear(
                grid.porosity(global), 0.25, 1.0e-12,
                "StructuredGrid porosity mismatch");
        }
    }

    Vec solution = grid.createGlobalVector(TestTag::numVars_);
    PetscInt vectorSize = 0;
    PetscCallAbort(PETSC_COMM_WORLD, VecGetSize(solution, &vectorSize));
    MPMC::Test::require(vectorSize == 24 * TestTag::numVars_,
                        "StructuredGrid global vector size mismatch");

    Mat matrix = grid.createMatrix(TestTag::numVars_);
    PetscInt rows = 0;
    PetscInt columns = 0;
    PetscCallAbort(PETSC_COMM_WORLD, MatGetSize(matrix, &rows, &columns));
    MPMC::Test::require(rows == vectorSize && columns == vectorSize,
                        "StructuredGrid matrix size mismatch");

    PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&matrix));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));

    PetscPrintf(PETSC_COMM_WORLD,
                "StructuredGrid integration test: ALL PASS\n");
}

void testSubcommunicatorReport()
{
    PetscMPIInt worldRank = 0;
    PetscCallMPIAbort(
        PETSC_COMM_WORLD, MPI_Comm_rank(PETSC_COMM_WORLD, &worldRank));

    MPI_Comm subcomm = MPI_COMM_NULL;
    const int color = worldRank % 2;
    PetscCallMPIAbort(
        PETSC_COMM_WORLD,
        MPI_Comm_split(PETSC_COMM_WORLD, color, worldRank, &subcomm));

    // 只有一个子通信域进入 collective 报告。若报告错误使用 WORLD，
    // 这里会与另一个颜色的 rank 不匹配并挂起。
    if (color == 0)
    {
        Grid grid(2, 1, 1, MPMC::CellSize{1.0, 1.0, 1.0}, subcomm);
        grid.setup();
        grid.initializeRockProperties();
        PetscCallAbort(
            subcomm, VecSet(grid.permeabilityVectorHandle(), 1.0e-13));
        PetscCallAbort(
            subcomm, VecSet(grid.porosityVectorHandle(), 0.2));
        MPMC::printStructuredGridSummary(grid);
    }

    PetscCallMPIAbort(
        PETSC_COMM_WORLD, MPI_Barrier(PETSC_COMM_WORLD));
    PetscCallMPIAbort(PETSC_COMM_WORLD, MPI_Comm_free(&subcomm));
}

} // namespace

int main(int argc, char **argv)
{
    const PetscErrorCode initError = PetscInitialize(&argc, &argv, nullptr, nullptr);
    if (initError != PETSC_SUCCESS)
    {
        return static_cast<int>(initError);
    }

    int status = 0;
    try
    {
        runTest();
        testSubcommunicatorReport();
    }
    catch (const std::exception &error)
    {
        PetscPrintf(PETSC_COMM_WORLD,
                    "StructuredGrid integration test failed: %s\n",
                    error.what());
        status = 1;
    }

    const PetscErrorCode finalError = PetscFinalize();
    return finalError == PETSC_SUCCESS ? status : static_cast<int>(finalError);
}
