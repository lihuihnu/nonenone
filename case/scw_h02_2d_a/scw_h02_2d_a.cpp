#include <case/petsc_custom_hooks.hpp>
#include <case/structured_single_eos_reservoir_runner.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <filesystem>
#include <fstream>

namespace Case
{
struct Definition
{
    using Config = H02A::Config;
    static auto wells() { return H02A::makeWellDefinitions(60, 20); }
};

using Runner = MPMC::cases::StructuredSingleEosReservoirRunner<Definition>;
using Indices = Runner::Indices;
using Runtime = Runner::Runtime;
using Grid = Runner::Grid;

static_assert(!Indices::fullyCompositionalThreePhase);
static_assert(Indices::hasIndependentWaterConservation);
static_assert(Indices::numComponents == 1);

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscInt nx = 60;
    PetscInt ny = 20;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetInt(nullptr, nullptr, "-h02_nx", &nx, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetInt(nullptr, nullptr, "-h02_ny", &ny, nullptr));
    if (nx < 2 || ny < 2)
        throw std::invalid_argument("H02 A requires nx,ny >= 2.");

    auto run = MPMC::cases::readRunOptions<H02A::Config>();
    MPMC::cases::validateCaseConfig<Indices, H02A::Config>();

    Grid grid(
        nx, ny, 1,
        MPMC::GridExtent{0.30, 0.10, 0.010});
    grid.setup();
    Runner::initializeRock(grid);

    auto fluid = MPMC::cases::makeFluidSystem<Indices, H02A::Config>();
    using Eval = typename Indices::ValueType;
    fluid.flowViscosityOverride =
        [](Eval, const MPMC::FluidSystem<Indices>::Composition &) {
            return Eval(2.0e-3);
        };

    auto opt =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, H02A::Config>(run);
    opt.scaling.pressureScale = 1.0e7;
    opt.scaling.compositionScale = 1.0;
    opt.scaling.saturationScale = 1.0;
    opt.scaling.massResidualScale = 2.0e-5;
    opt.scaling.fugacityResidualScale = 1.0;
    opt.scaling.closureResidualScale = 1.0;
    opt.scaling.rateWellResidualFloor = 2.5e-8;
    MPMC::cases::applyNaturalScalingPetscOptions(opt, rank);

    Runtime runtime(grid, fluid, opt);
    auto definitions = H02A::makeWellDefinitions(nx, ny);
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid,
        definitions,
        [](Grid &g, const auto &d, int i, int j, int k) {
            return Runner::wellIndex(g, i, j, k, d.radius, d.skin);
        });
    runtime.setWells(std::move(wells));

    Vec x = grid.createGlobalVector(Indices::numPrimaryVariables);
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(x, 0.0));
    auto values = grid.vecGetArray<Indices::numPrimaryVariables>(x);
    const auto owned = grid.ownedRegion();
    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                auto &p = values[k][j][i];
                p.fill(0.0);
                p[Indices::Primary::pressure] = 28.0e6;
                p[Indices::Primary::waterSaturation] = 0.0;
                p[Indices::Primary::liquidSaturation] = 1.0;
                p[Indices::Primary::vaporSaturation] = 0.0;
            }

    for (const auto &d : definitions)
    {
        const int i = d.completion.i;
        const int j = d.completion.j;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            owned.zStart == 0)
        {
            values[0][j][i][Indices::Primary::wellPressure] = d.initialBhp;
        }
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(x, values);

    runtime.initializePhaseStateFromSolution(
        x, MPMC::HydrocarbonPhaseState::LiquidOnly);
    runtime.updateState(x);
    runtime.initializeHistory(x);

    if (rank == 0)
    {
        std::filesystem::create_directories(run.resultDirectory);
        std::ofstream manifest(run.resultDirectory + "/h02_manifest.json");
        manifest
            << "{\"mode\":\"A\",\"nx\":" << nx
            << ",\"ny\":" << ny
            << ",\"PV_m3\":7.5e-5"
            << ",\"Q_target_m3_s\":2.5e-8"
            << ",\"physics\":\"IMMISCIBLE_INDEPENDENT_WATER_PLUS_PURE_HEAVY_PR\""
            << ",\"target_pvi\":" << run.targetPVI
            << "}\n";
    }

    MPMC::cases::NaturalSolver<Runtime> solver(
        runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, H02A::Config>(
        runtime, solver.snes(), x, run, rank);

    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
    return 0;
}
}

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return MPMC::cases::runPetscCaseMain(
        argc, argv, [] { return Case::run(); });
}
