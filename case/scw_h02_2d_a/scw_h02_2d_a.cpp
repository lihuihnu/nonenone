#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/case_support.hpp>
#include <case/natural_scaling_options.hpp>
#include <case/petsc_case_main.hpp>
#include <case/well_factory.hpp>
#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/phase_state.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace Case
{
using Config = H02A::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    H02A::Model::numberOfComponents,
    H02A::Model::hasWater,
    H02A::Model::hasWells,
    H02A::Model::enableDissolution,
    H02A::Model::enableAdsorption,
    H02A::Model::enableLandTrapping,
    H02A::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;

static_assert(!Indices::fullyCompositionalThreePhase);
static_assert(Indices::hasIndependentWaterConservation);
static_assert(Indices::numComponents == 1);

void initializeRock(Grid &grid)
{
    grid.initializeRockProperties();
    auto permeability = grid.vecGetArray<3>(grid.permeabilityVector);
    auto porosity = grid.vecGetArray<1>(grid.porosityVector);
    const auto owned = grid.ownedRegion();
    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                permeability[k][j][i] = {
                    H02A::Rock::kx, H02A::Rock::ky, H02A::Rock::kz};
                porosity[k][j][i][0] = H02A::Rock::porosity;
            }
    grid.vecRestoreArray<3>(grid.permeabilityVector, permeability);
    grid.vecRestoreArray<1>(grid.porosityVector, porosity);
}

MPMC::FluidSystem<Indices> makeFluid()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, Config>();
    using Eval = typename Indices::ValueType;

    // A is a true no-transfer counterfactual: Water is an independent
    // conserved phase.  The Heavy EOS remains pure PR; its flow viscosity is
    // the same 2 mPa.s reference used by B before composition feedback.
    fluid.flowViscosityOverride =
        [](Eval, const MPMC::FluidSystem<Indices>::Composition &) {
            return Eval(2.0e-3);
        };
    fluid.waterViscosity = [](Eval) { return Eval(5.0e-5); };
    fluid.waterFormationVolumeFactor = [](Eval) { return Eval(1.0); };
    return fluid;
}

std::vector<Well> makeWells(Grid &grid, int nx, int ny)
{
    const auto definitions = H02A::makeWellDefinitions(nx, ny);
    return MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid,
        definitions,
        [](Grid &g, const auto &definition, int i, int j, int k) {
            const auto size = g.cellSize({i,j,k});
            return MPMC::verticalPeacemanWellIndex({
                size[0], size[1], size[2],
                H02A::Rock::kx, H02A::Rock::ky,
                definition.radius, definition.skin});
        });
}

void initializeSolution(
    Grid &grid,
    Vec solution,
    int nx,
    int ny)
{
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(solution, 0.0));
    auto values = grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();

    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                auto &primary = values[k][j][i];
                primary.fill(0.0);
                primary[Indices::Primary::pressure] = 28.0e6;
                primary[Indices::Primary::liquidSaturation] = 1.0;
                primary[Indices::Primary::vaporSaturation] = 0.0;
                primary[Indices::Primary::waterSaturation] = 0.0;
            }

    for (const auto &definition : H02A::makeWellDefinitions(nx, ny))
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            owned.zStart <= 0 && 0 < owned.zStart + owned.zCount)
        {
            values[0][j][i][Indices::Primary::wellPressure] =
                definition.initialBhp;
        }
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscInt nx = H02A::Grid::nx;
    PetscInt ny = H02A::Grid::ny;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetInt(nullptr, nullptr, "-h02_nx", &nx, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetInt(nullptr, nullptr, "-h02_ny", &ny, nullptr));
    if (nx < 2 || ny < 2)
        throw std::invalid_argument("H02 A requires nx,ny >= 2.");

    MPMC::cases::validateCaseConfig<Indices, Config>();
    auto run = MPMC::cases::readRunOptions<Config>();
    MPMC::cases::printRunSummary(Config::name, run, rank);

    Grid grid(
        nx, ny, 1,
        MPMC::GridExtent{H02A::Grid::lx, H02A::Grid::ly, H02A::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = makeFluid();
    auto runtimeOptions =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);
    runtimeOptions.scaling.pressureScale = 1.0e7;
    runtimeOptions.scaling.compositionScale = 1.0;
    runtimeOptions.scaling.saturationScale = 1.0;
    runtimeOptions.scaling.massResidualScale = 2.0e-5;
    runtimeOptions.scaling.fugacityResidualScale = 1.0;
    runtimeOptions.scaling.closureResidualScale = 1.0;
    runtimeOptions.scaling.rateWellResidualFloor = 2.5e-8;
    MPMC::cases::applyNaturalScalingPetscOptions(runtimeOptions, rank);

    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid, nx, ny));

    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, solution, nx, ny);
    runtime.initializePhaseStateFromSolution(
        solution, MPMC::HydrocarbonPhaseState::LiquidOnly);
    runtime.updateState(solution);
    runtime.initializeHistory(solution);

    if (rank == 0)
    {
        std::filesystem::create_directories(run.resultDirectory);
        std::ofstream manifest(run.resultDirectory + "/h02_manifest.json");
        manifest
            << "{\"mode\":\"A\",\"nx\":" << nx
            << ",\"ny\":" << ny
            << ",\"PV_m3\":7.5e-5"
            << ",\"Q_target_m3_s\":2.5e-8"
            << ",\"physics\":\"INDEPENDENT_WATER_NO_COMPONENT_TRANSFER\""
            << ",\"target_pvi\":" << run.targetPVI
            << "}\n";
    }

    MPMC::cases::NaturalSolver<Runtime> solver(
        runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, Config>(
        runtime, solver.snes(), solution, run, rank);

    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));
    return 0;
}
} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return MPMC::cases::runPetscCaseMain(
        argc, argv, [] { return Case::run(); });
}
