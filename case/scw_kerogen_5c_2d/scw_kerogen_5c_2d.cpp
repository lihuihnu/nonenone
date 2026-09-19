#include <case/petsc_custom_hooks.hpp>
#include <case/structured_single_eos_reservoir_runner.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <filesystem>
#include <fstream>

namespace Case {

struct Definition {
    using Config = ScwKerogen5C2D::Config;
    static auto wells() { return ScwKerogen5C2D::makeWellDefinitions(60, 20); }
};

using Runner = MPMC::cases::StructuredSingleEosReservoirRunner<Definition>;
using Indices = Runner::Indices;
using Runtime = Runner::Runtime;
using Grid = Runner::Grid;

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscInt nx = 60, ny = 20;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetInt(nullptr, nullptr, "-mc_nx", &nx, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetInt(nullptr, nullptr, "-mc_ny", &ny, nullptr));
    if (nx < 2 || ny < 2)
        throw std::invalid_argument("5C screening requires nx,ny >= 2");

    auto run = MPMC::cases::readRunOptions<ScwKerogen5C2D::Config>();
    MPMC::cases::validateCaseConfig<Indices, ScwKerogen5C2D::Config>();

    Grid grid(
        nx, ny, 1,
        MPMC::GridExtent{0.30, 0.10, 0.010});
    grid.setup();
    Runner::initializeRock(grid);

    auto fluid =
        MPMC::cases::makeFluidSystem<Indices, ScwKerogen5C2D::Config>();

    auto opt =
        MPMC::cases::makeRuntimeOptions<
            Indices, Runtime, ScwKerogen5C2D::Config>(run);
    opt.scaling.pressureScale = 1e7;
    opt.scaling.compositionScale = 1.0;
    opt.scaling.saturationScale = 1.0;
    opt.scaling.massResidualScale = 2e-5;
    opt.scaling.fugacityResidualScale = 1.0;
    opt.scaling.closureResidualScale = 1.0;
    opt.scaling.rateWellResidualFloor = 2.5e-8;
    MPMC::cases::applyNaturalScalingPetscOptions(opt, rank);

    Runtime runtime(grid, fluid, opt);
    const auto definitions =
        ScwKerogen5C2D::makeWellDefinitions(nx, ny);
    auto wells =
        MPMC::cases::makeStructuredWells<Indices, PetscInt>(
            grid, definitions,
            [](Grid &g, const auto &d, int i, int j, int k) {
                return Runner::wellIndex(
                    g, i, j, k, d.radius, d.skin);
            });
    runtime.setWells(std::move(wells));

    Vec x = grid.createGlobalVector(Indices::numPrimaryVariables);
    runtime.initializeUniformFromPTZ(
        x,
        ScwKerogen5C2D::InitialState::pressure,
        ScwKerogen5C2D::InitialState::temperature,
        ScwKerogen5C2D::InitialState::overallComposition);

    auto array =
        grid.vecGetArray<Indices::numPrimaryVariables>(x);
    const auto own = grid.ownedRegion();
    for (const auto &definition : definitions)
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        if (i >= own.xStart && i < own.xStart + own.xCount &&
            j >= own.yStart && j < own.yStart + own.yCount &&
            own.zStart == 0)
        {
            array[0][j][i][Indices::Primary::wellPressure] =
                definition.initialBhp;
        }
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(x, array);
    runtime.initializeHistory(x);

    if (rank == 0)
    {
        std::filesystem::create_directories(run.resultDirectory);
        std::ofstream manifest(
            run.resultDirectory + "/multicomponent_manifest.json");
        manifest
            << "{\n"
            << "  \"scope\": \"CONDITIONAL_5C_FLOW_SCREENING_NOT_FORMAL_VALIDATION\",\n"
            << "  \"temperature_K\": 653.15,\n"
            << "  \"pressure_center_MPa\": 28.0,\n"
            << "  \"nx\": " << nx << ",\n"
            << "  \"ny\": " << ny << ",\n"
            << "  \"components\": [\"H2O\",\"OIL_GASOLINE\",\"OIL_DIESEL\",\"OIL_MIDDLE\",\"OIL_HEAVY\"],\n"
            << "  \"initial_z_H2O\": 0.20,\n"
            << "  \"target_pvi\": " << run.targetPVI << ",\n"
            << "  \"PV_m3\": 7.5e-5,\n"
            << "  \"Q_target_m3_s\": 2.5e-8,\n"
            << "  \"eos\": \"PR76\",\n"
            << "  \"bip_status\": \"SCREENING_PRIORS_NOT_FINAL_CALIBRATION\"\n"
            << "}\n";
    }

    MPMC::cases::NaturalSolver<Runtime> solver(
        runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<
        Indices, Runtime, ScwKerogen5C2D::Config>(
            runtime, solver.snes(), x, run, rank);

    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
    return 0;
}

} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return MPMC::cases::runPetscCaseMain(
        argc, argv, [] { return Case::run(); });
}
