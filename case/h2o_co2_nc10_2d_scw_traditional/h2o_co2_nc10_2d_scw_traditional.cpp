/**
 * @file h2o_co2_nc10_2d_scw_traditional.cpp
 * @brief 653.15 K、28 MPa 超临界水 + PR(CO2-nC10) 算例。
 */
#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/case_support.hpp>
#include <case/well_factory.hpp>
#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/phase_state.hpp>
#include <natural/thermo/aqueous_volume.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <cmath>
#include <exception>
#include <vector>

namespace Case
{
using Config = TraditionalConfig::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    TraditionalConfig::Model::numberOfComponents, TraditionalConfig::Model::hasWater,
    TraditionalConfig::Model::hasWells, TraditionalConfig::Model::enableDissolution,
    TraditionalConfig::Model::enableAdsorption, TraditionalConfig::Model::enableLandTrapping,
    TraditionalConfig::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;

static_assert(!Indices::fullyCompositionalThreePhase);
static_assert(Indices::hasIndependentWaterConservation);
static_assert(BenchmarkCommon::temperature >
              BenchmarkCommon::waterCriticalTemperature);
static_assert(BenchmarkCommon::initialPressure >
              BenchmarkCommon::waterCriticalPressure);
static_assert(BenchmarkCommon::baseWells[0].maximumBhp ==
              BenchmarkCommon::injectorMaximumBhp);

MPMC::FluidSystem<Indices> makeFluid()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, TraditionalConfig::PrFactoryConfig>();
    BenchmarkCommon::applyCommonRelativePermeability(fluid);
    using Eval = typename Indices::ValueType;
    constexpr double referenceDensity = TraditionalConfig::Fluid::surfaceDensity[2];
    fluid.waterFormationVolumeFactor = [=](Eval pressure) {
        const Eval density = MPMC::IapwsIf97WaterDensity::density(
            pressure, Eval(BenchmarkCommon::temperature));
        return density / Eval(referenceDensity);
    };
    fluid.waterViscosity = [](Eval pressure) {
        constexpr double a = -3.705013;
        constexpr double b = 0.00289258;
        constexpr double c = 3.98950;
        constexpr double d = -0.00326;
        constexpr double t0 = 141.5;
        const Eval reducedPressure = pressure / 1.0e6;
        const double shiftedTemperature = BenchmarkCommon::temperature / t0 - 1.0;
        return exp(a + b * reducedPressure + (c + d * reducedPressure) / shiftedTemperature)
            * 1.0e-3;
    };
    return fluid;
}

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
                    TraditionalConfig::Rock::kx(i,j,k), TraditionalConfig::Rock::ky(i,j,k),
                    TraditionalConfig::Rock::kz(i,j,k)};
                porosity[k][j][i][0] = TraditionalConfig::Rock::porosity(i,j,k);
            }
    grid.vecRestoreArray<3>(grid.permeabilityVector, permeability);
    grid.vecRestoreArray<1>(grid.porosityVector, porosity);
}

std::vector<Well> makeWells(Grid &grid)
{
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid, TraditionalWellConfig::wells,
        [](Grid &g, const auto &definition, int i, int j, int k) {
            const auto size = g.cellSize({i,j,k});
            return MPMC::verticalPeacemanWellIndex({
                size[0], size[1], size[2], TraditionalConfig::Rock::kx(i,j,k),
                TraditionalConfig::Rock::ky(i,j,k), definition.radius, definition.skin});
        });
    wells[0].control = MPMC::WellControl::ReservoirTotalRate;
    wells[0].primaryControl = MPMC::WellControl::ReservoirTotalRate;
    return wells;
}

void initializeSolution(Grid &grid, Vec solution)
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
                primary[Indices::Primary::pressure] = BenchmarkCommon::initialPressure;
                primary[Indices::Primary::liquidComposition[0]] = 0.0;
                primary[Indices::Primary::vaporComposition[0]] = 0.0;
                primary[Indices::Primary::liquidSaturation] = BenchmarkCommon::targetOilSaturation;
                primary[Indices::Primary::vaporSaturation] = 0.0;
                primary[Indices::Primary::waterSaturation] = BenchmarkCommon::targetWaterSaturation;
            }
    for (const auto &definition : TraditionalWellConfig::wells)
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        const int k = definition.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
            values[k][j][i][Indices::Primary::wellPressure] = definition.initialBhp;
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

void configureBenchmarkSolverDefaults()
{
    PetscBool lineSearchSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-snes_linesearch_type", &lineSearchSet));
    if (!lineSearchSet)
        PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(
            nullptr, "-snes_linesearch_type", "basic"));

    PetscBool stepToleranceSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-snes_stol", &stepToleranceSet));
    if (!stepToleranceSet)
        PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(
            nullptr, "-snes_stol", "1e-12"));
}

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPMC::cases::validateCaseConfig<Indices, Config>();
    auto run = MPMC::cases::readRunOptions<Config>();
    MPMC::cases::printRunSummary(TraditionalConfig::name, run, rank);
    Grid grid(TraditionalConfig::Grid::nx, TraditionalConfig::Grid::ny,
              TraditionalConfig::Grid::nz,
              MPMC::GridExtent{TraditionalConfig::Grid::lx,
                               TraditionalConfig::Grid::ly,
                               TraditionalConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);
    auto fluid = makeFluid();
    if (rank == 0)
        PetscPrintf(PETSC_COMM_SELF,
            "[BENCHMARK] model=traditional grid=60x20x1 T=%.8g K P=%.8g Pa "
            "PV=%.12g m3 rate=%.12g m3/s "
            "inj_max_bhp=%.12g Pa prod_bhp=%.12g Pa\n"
            "[INIT] independent-H2O + PR(CO2/nC10), S(O/G/W)=0.8/0/0.2\n",
            BenchmarkCommon::temperature, BenchmarkCommon::initialPressure,
            BenchmarkCommon::poreVolume, BenchmarkCommon::injectionRate,
            BenchmarkCommon::injectorMaximumBhp, BenchmarkCommon::producerBhp);

    auto runtimeOptions = MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);
    runtimeOptions.scaling.enabled = true;
    runtimeOptions.scaling.rateWellResidualFloor = BenchmarkCommon::injectionRate;
    runtimeOptions.scaling.validate();
    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid));
    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, solution);
    runtime.initializePhaseStateFromSolution(solution, MPMC::HydrocarbonPhaseState::LiquidOnly);
    runtime.updateState(solution);
    runtime.initializeHistory(solution);
    configureBenchmarkSolverDefaults();
    MPMC::cases::NaturalSolver<Runtime> solver(runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, Config>(
        runtime, solver.snes(), solution, run, rank);
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));
    return 0;
}
} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    PetscCallAbort(PETSC_COMM_WORLD, PetscInitialize(&argc, &argv, nullptr, nullptr));
    int code = 0;
    try { code = Case::run(); }
    catch (const std::exception &error) {
        PetscPrintf(PETSC_COMM_WORLD, "[ERROR][CASE] %s\n", error.what());
        code = 1;
    }
    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
    return code;
}
