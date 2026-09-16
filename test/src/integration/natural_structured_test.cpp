/**
 * @file natural_structured_test.cpp
 * @brief 集成测试：验证 `natural_structured_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/petsc/callbacks.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace
{

using Config =
    MPMC::CompositionalModelConfig<
        6,
        true,   // water
        false,  // no well unknown in this minimal runtime smoke test
        false,
        false,
        false>;

using Indices = MPMC::ADIndices<Config>;
using Grid = MPMC::StructuredGridCore;
using Runtime =
    MPMC::NaturalStructuredGridRuntime<
        Indices,
        Grid>;

constexpr double mdToM2 = 9.869233e-16;

MPMC::CompositionalMixture<Indices> makeMixture()
{
    return MPMC::CompositionalMixture<Indices>(
        {189.515, 304.2, 387.607, 597.497, 698.515, 875.0},
        {4580011.59, 7386592.50, 4095515.97,
         3345244.875, 1768374.5625, 1169006.79},
        {9.97012032965401e-5, 9.26344713533338e-5,
         0.000217076707259486, 0.000381162235869935,
         0.000721410148917871, 0.00113570073874421},
        {0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141},
        {0.0161594, 0.04401, 0.0455725, 0.11774, 0.248827, 0.48152},
        {{0.0, 0.00070981, 0.00077754, 0.0100, 0.0110, 0.0110},
         {0.00070981, 0.0, 0.1500, 0.1500, 0.1500, 0.1500},
         {0.00077754, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0100, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0110, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0110, 0.1500, 0.0, 0.0, 0.0, 0.0}});
}

MPMC::FluidSystem<Indices> makeFluidSystem()
{
    MPMC::CubicEquationOfState<Indices> eos(
        0.4572355,
        0.0779691,
        makeMixture(),
        1,
        2.414213562373095,
        -0.414213562373095);

    MPMC::FluidSystem<Indices> fluid(
        {800.0, 2.0, 1000.0},
        {1.0e-3, 1.0e-5, 2.0e-4},
        std::move(eos),
        {"N2/CH4", "CO2", "C2-5", "C6-13", "C14-24", "C25-80"},
        {"liquid", "vapor", "water"},
        387.45);

    using Eval = typename Indices::ValueType;

    fluid.gasRelativePermeability = [](Eval s) { return s * s; };
    fluid.waterRelativePermeability = [](Eval s) { return s * s; };
    fluid.threePhaseOilRelativePermeability =
        [](Eval, Eval so, Eval) { return so * so; };
    fluid.waterViscosity = [](Eval) { return Eval(2.0e-4); };
    fluid.waterFormationVolumeFactor = [](Eval) { return Eval(1.0); };

    return fluid;
}

void initializeRock(Grid &grid)
{
    grid.initializeRockProperties();

    auto rock = grid.rockWriteView();
    const auto owned = grid.ownedRegion();

    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                rock.permeability(i, j, k) = {
                    50.0 * mdToM2,
                    50.0 * mdToM2,
                    50.0 * mdToM2};
                rock.porosity(i, j, k) = 0.25;
            }
}

void fillSolution(Grid &grid, Vec solution)
{
    auto values =
        grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();

    const std::array<double, 5> composition{
        0.3246914,
        0.0128351,
        0.2278401,
        0.2606985,
        0.1134144};

    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                auto &x = values[k][j][i];
                x.fill(0.0);

                const int cell = grid.ijkToGlobal(i, j, k);
                x[Indices::Primary::pressure] =
                    15.0e6 - 100.0 * static_cast<double>(cell);

                for (int c = 0;
                     c < Indices::numIndependentCompositionsPerPhase;
                     ++c)
                {
                    x[Indices::Primary::liquidComposition[
                        static_cast<std::size_t>(c)]] =
                        composition[static_cast<std::size_t>(c)];
                    x[Indices::Primary::vaporComposition[
                        static_cast<std::size_t>(c)]] =
                        composition[static_cast<std::size_t>(c)];
                }

                x[Indices::Primary::waterSaturation] = 0.20;
                x[Indices::Primary::liquidSaturation] = 0.50;
                x[Indices::Primary::vaporSaturation] = 0.30;
            }

    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

void runTest()
{
    Grid grid(4, 3, 2, MPMC::GridExtent{40.0, 30.0, 10.0});
    grid.setup();
    initializeRock(grid);

    auto fluid = makeFluidSystem();
    Runtime runtime(grid, fluid);

    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    fillSolution(grid, solution);

    runtime.initializePhaseStateFromSolution(solution);
    runtime.updateState(solution);
    runtime.initializeHistory(solution);

    Vec residual = runtime.createResidualVector();
    Mat jacobian = runtime.createJacobian();

    SNES snes = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        MPMC::installNaturalCallbacks(snes, runtime, residual, jacobian));
    PetscCallAbort(PETSC_COMM_WORLD, SNESSetFromOptions(snes));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        MPMC::installNaturalLineSearchCallbacks(snes, runtime));

    SNESLineSearch lineSearch = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, SNESGetLineSearch(snes, &lineSearch));
    PetscErrorCode (*preCheck)(
        SNESLineSearch, Vec, Vec, PetscBool *, void *) = nullptr;
    PetscErrorCode (*postCheck)(
        SNESLineSearch, Vec, Vec, Vec, PetscBool *, PetscBool *, void *) = nullptr;
    void *preContext = nullptr;
    void *postContext = nullptr;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        SNESLineSearchGetPreCheck(lineSearch, &preCheck, &preContext));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        SNESLineSearchGetPostCheck(lineSearch, &postCheck, &postContext));
    PetscCheckAbort(
        preCheck != nullptr && postCheck != nullptr &&
            preContext == &runtime && postContext == &runtime,
        PETSC_COMM_WORLD,
        PETSC_ERR_PLIB,
        "Natural standard PETSc line-search callbacks are not installed.");

    Vec step = nullptr;
    Vec candidate = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(solution, &step));
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(solution, &candidate));
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(step, 0.0));
    PetscCallAbort(PETSC_COMM_WORLD, VecCopy(solution, candidate));
    PetscBool changedStep = PETSC_FALSE;
    PetscBool changedDirection = PETSC_FALSE;
    PetscBool changedCandidate = PETSC_FALSE;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        SNESLineSearchPreCheck(
            lineSearch, solution, step, &changedStep));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        SNESLineSearchPostCheck(
            lineSearch, solution, step, candidate,
            &changedDirection, &changedCandidate));
    PetscCheckAbort(
        changedStep == PETSC_TRUE && changedCandidate == PETSC_TRUE,
        PETSC_COMM_WORLD,
        PETSC_ERR_PLIB,
        "Natural standard PETSc line-search callbacks were not executed.");

    PetscCallAbort(
        PETSC_COMM_WORLD,
        runtime.formFunction(solution, residual));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        runtime.formJacobian(solution, jacobian, jacobian));

    PetscReal residualNorm = 0.0;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        VecNorm(residual, NORM_2, &residualNorm));

    PetscCheckAbort(
        std::isfinite(static_cast<double>(residualNorm)),
        PETSC_COMM_WORLD,
        PETSC_ERR_FP,
        "Natural StructuredGrid residual norm is not finite.");

    const auto diagnostics = runtime.evaluateGlobalDiagnostics(solution);
    PetscCheckAbort(
        diagnostics.cellCount == 24,
        PETSC_COMM_WORLD,
        PETSC_ERR_PLIB,
        "Natural StructuredGrid diagnostics cell count mismatch.");
    PetscCheckAbort(
        std::abs(diagnostics.porosityAverage - 0.25) < 1.0e-12 &&
        std::abs(diagnostics.saturationAverage[Indices::Phase::liquid] - 0.50) < 1.0e-12 &&
        std::abs(diagnostics.saturationAverage[Indices::Phase::vapor] - 0.30) < 1.0e-12 &&
        std::abs(diagnostics.saturationAverage[Indices::Phase::water] - 0.20) < 1.0e-12,
        PETSC_COMM_WORLD,
        PETSC_ERR_PLIB,
        "Natural StructuredGrid diagnostics averages mismatch.");

    PetscInt m = 0, n = 0;
    PetscCallAbort(PETSC_COMM_WORLD, MatGetSize(jacobian, &m, &n));

    const PetscInt expected =
        4 * 3 * 2 * Indices::numPrimaryVariables;

    PetscCheckAbort(
        m == expected && n == expected,
        PETSC_COMM_WORLD,
        PETSC_ERR_PLIB,
        "Natural StructuredGrid Jacobian size mismatch.");

    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscPrintf(
            PETSC_COMM_WORLD,
            "Natural/StructuredGrid runtime test: ALL PASS, residual norm = %.6e\n",
            static_cast<double>(residualNorm)));

    PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&jacobian));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&residual));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&candidate));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&step));
    PetscCallAbort(PETSC_COMM_WORLD, SNESDestroy(&snes));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));
}

} // namespace

int main(int argc, char **argv)
{
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscInitialize(&argc, &argv, nullptr, nullptr));

    try
    {
        runTest();
    }
    catch (const std::exception &error)
    {
        PetscPrintf(PETSC_COMM_WORLD, "FAILED: %s\n", error.what());
        PetscFinalize();
        return 1;
    }

    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
    return 0;
}
