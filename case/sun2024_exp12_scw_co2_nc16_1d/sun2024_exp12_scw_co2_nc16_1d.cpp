/**
 * @file sun2024_exp12_scw_co2_nc16_1d.cpp
 * @brief Sun-2024 Exp.12 派生独立超临界水/CO2/nC16 一维驱替入口。
 */
#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "case_fluid.hpp"
#include "well_config.hpp"

#include <case/case_support.hpp>
#include <case/well_factory.hpp>
#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/thermo/aqueous_volume.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace Case
{
using Config = CaseConfig::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents, CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells, CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption, CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;

static_assert(!Indices::fullyCompositionalThreePhase);
static_assert(Indices::hasIndependentWaterConservation);
static_assert(CaseConfig::Grid::nx * CaseConfig::Grid::ny *
                  CaseConfig::Grid::nz == 48);
static_assert(Sun2024Exp12::temperature > 647.096);
static_assert(Sun2024Exp12::initialPressure > 22.064e6);

MPMC::FluidSystem<Indices> makeFluid()
{
    auto fluid = MPMC::cases::makeFluidSystem<
        Indices, CaseConfig::PrFactoryConfig>();
    BenchmarkCaseFluid::apply(fluid);

    using Eval = typename Indices::ValueType;
    fluid.waterFormationVolumeFactor = [](Eval pressure) {
        const Eval reservoirDensity = MPMC::IapwsIf97WaterDensity::density(
            pressure, Eval(Sun2024Exp12::temperature));
        return reservoirDensity /
            Eval(CaseConfig::Fluid::surfaceDensity[Indices::Phase::water]);
    };
    fluid.waterViscosity = [](Eval pressure) {
        constexpr double a = -3.705013;
        constexpr double b = 0.00289258;
        constexpr double c = 3.98950;
        constexpr double d = -0.00326;
        constexpr double t0 = 141.5;
        const Eval reducedPressure = pressure / 1.0e6;
        const double shiftedTemperature = Sun2024Exp12::temperature / t0 - 1.0;
        return exp(a + b * reducedPressure +
                   (c + d * reducedPressure) / shiftedTemperature) * 1.0e-3;
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
                    CaseConfig::Rock::kx(i, j, k),
                    CaseConfig::Rock::ky(i, j, k),
                    CaseConfig::Rock::kz(i, j, k)};
                porosity[k][j][i][0] = CaseConfig::Rock::porosity(i, j, k);
            }
    grid.vecRestoreArray<3>(grid.permeabilityVector, permeability);
    grid.vecRestoreArray<1>(grid.porosityVector, porosity);
}

std::vector<Well> makeWells(Grid &grid)
{
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid, WellConfig::wells,
        [](Grid &g, const auto &definition, int i, int j, int k) {
            const auto size = g.cellSize({i, j, k});
            return MPMC::verticalPeacemanWellIndex({
                size[0], size[1], size[2], CaseConfig::Rock::kx(i, j, k),
                CaseConfig::Rock::ky(i, j, k), definition.radius,
                definition.skin});
        });

    auto &injector = wells[0];
    injector.control = MPMC::WellControl::ReservoirTotalRate;
    injector.primaryControl = MPMC::WellControl::ReservoirTotalRate;
    injector.injectionPhaseFraction.fill(0.0);
    injector.injectionPhaseFraction[Indices::Phase::water] =
        Sun2024Exp12::waterInjectionVolumeFraction;
    injector.injectionPhaseFraction[Indices::Phase::vapor] =
        Sun2024Exp12::co2InjectionVolumeFraction;
    injector.injectionComponentMassFraction.fill(0.0);
    injector.injectionComponentMassFraction[CaseConfig::Fluid::co2Component] =
        1.0;
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
                primary[Indices::Primary::pressure] =
                    CaseConfig::InitialState::pressure;
                primary[Indices::Primary::liquidComposition[0]] =
                    CaseConfig::InitialState::oilComposition[0];
                primary[Indices::Primary::vaporComposition[0]] =
                    CaseConfig::InitialState::gasComposition[0];
                primary[Indices::Primary::liquidSaturation] =
                    CaseConfig::InitialState::oilSaturation;
                primary[Indices::Primary::vaporSaturation] =
                    CaseConfig::InitialState::gasSaturation;
                primary[Indices::Primary::waterSaturation] =
                    CaseConfig::InitialState::waterSaturation;
            }

    for (const auto &definition : WellConfig::wells)
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        const int k = definition.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
            values[k][j][i][Indices::Primary::wellPressure] =
                definition.initialBhp;
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

void configureSolverDefaults()
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
    char requested[PETSC_MAX_PATH_LEN]{};
    PetscBool resultSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetString(
        nullptr, nullptr, "-result_dir", requested, sizeof(requested),
        &resultSet));
    if (!resultSet)
        run.resultDirectory = "./results/new-pr";
    MPMC::cases::printRunSummary(CaseConfig::name, run, rank);

    Grid grid(CaseConfig::Grid::nx, CaseConfig::Grid::ny,
              CaseConfig::Grid::nz,
              MPMC::GridExtent{CaseConfig::Grid::lx, CaseConfig::Grid::ly,
                               CaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = makeFluid();
    const double waterDensity = MPMC::IapwsIf97WaterDensity::density(
        Sun2024Exp12::initialPressure, Sun2024Exp12::temperature);
    if (rank == 0)
        PetscPrintf(PETSC_COMM_SELF,
            "[BENCHMARK] source=Sun2024-Exp12-derived model=new-pr "
            "grid=48x1x1 T=%.8g K P=%.8g Pa PV=%.12g m3 horizon=%.8g PV\n"
            "[SCW] independent mobile phase, IF97 density=%.12g kg/m3\n"
            "[INJECTION] q(H2O/CO2)=%.12g/%.12g m3/s; total=%.12g m3/s\n"
            "[INIT] pure-nC16 oil + SCW, S(O/G/W)=%.12g/0/%.12g\n",
            Sun2024Exp12::temperature, Sun2024Exp12::initialPressure,
            Sun2024Exp12::poreVolume,
            Sun2024Exp12::targetInjectedPoreVolumes, waterDensity,
            Sun2024Exp12::waterInjectionRate,
            Sun2024Exp12::co2InjectionRate,
            Sun2024Exp12::totalInjectionRate,
            Sun2024Exp12::targetOilSaturation,
            Sun2024Exp12::targetSupercriticalFluidSaturation);

    auto runtimeOptions =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);
    runtimeOptions.scaling.enabled = true;
    runtimeOptions.scaling.rateWellResidualFloor =
        Sun2024Exp12::totalInjectionRate;
    runtimeOptions.scaling.validate();

    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid));
    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, solution);
    runtime.initializePhaseStateFromSolution(
        solution, MPMC::HydrocarbonPhaseState::LiquidOnly);
    runtime.updateState(solution);
    runtime.initializeHistory(solution);
    configureSolverDefaults();
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
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscInitialize(&argc, &argv, nullptr, nullptr));
    int code = 0;
    try
    {
        code = Case::run();
    }
    catch (const std::exception &error)
    {
        PetscPrintf(PETSC_COMM_WORLD, "[ERROR][CASE] %s\n", error.what());
        code = 1;
    }
    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
    return code;
}
