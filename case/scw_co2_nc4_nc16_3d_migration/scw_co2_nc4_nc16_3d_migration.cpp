/**
 * @file scw_co2_nc4_nc16_3d_migration.cpp
 * @brief 三维超临界水/CO2/nC4/nC16 运移与对照试验入口。
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

#include <cmath>
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
                  CaseConfig::Grid::nz == 2304);
static_assert(ScwMigration3D::temperature >
                  ScwMigration3D::waterCriticalTemperature);
static_assert(ScwMigration3D::initialPressure >
                  ScwMigration3D::waterCriticalPressure);

double readScwInjectionFraction()
{
    PetscReal value = ScwMigration3D::defaultScwInjectionVolumeFraction;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-scw_fraction", &value, nullptr));
    if (!std::isfinite(value) || value < 0.0 || value > 1.0)
        throw std::invalid_argument("-scw_fraction must be in [0,1].");
    return value;
}

MPMC::FluidSystem<Indices> makeFluid()
{
    auto fluid = MPMC::cases::makeFluidSystem<
        Indices, CaseConfig::PrFactoryConfig>();
    ScwMigrationCaseFluid::apply(fluid);
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
                    CaseConfig::Rock::kx(i, j, k),
                    CaseConfig::Rock::ky(i, j, k),
                    CaseConfig::Rock::kz(i, j, k)};
                rock.porosity(i, j, k) =
                    CaseConfig::Rock::porosity(i, j, k);
            }
}

std::vector<Well> makeWells(Grid &grid, double scwFraction)
{
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid, WellConfig::wells,
        [](Grid &wellGrid, const auto &definition, int i, int j, int k) {
            const auto size = wellGrid.cellSize({i, j, k});
            return MPMC::verticalPeacemanWellIndex({
                size[0], size[1], size[2],
                CaseConfig::Rock::kx(i, j, k),
                CaseConfig::Rock::ky(i, j, k),
                definition.radius, definition.skin});
        });

    auto &injector = wells[0];
    injector.control = MPMC::WellControl::ReservoirTotalRate;
    injector.primaryControl = MPMC::WellControl::ReservoirTotalRate;
    injector.injectionPhaseFraction.fill(0.0);
    injector.injectionPhaseFraction[Indices::Phase::water] = scwFraction;
    injector.injectionPhaseFraction[Indices::Phase::vapor] = 1.0 - scwFraction;
    injector.injectionComponentMassFraction.fill(0.0);
    injector.injectionComponentMassFraction[
        CaseConfig::Fluid::co2Component] = 1.0;

    // Balance the closed model with the same in-situ total volume rate at the
    // producer.  This avoids phase crossflow caused by a BHP sink/source while
    // preserving an auditable equal-throughput SCW-versus-CO2 comparison.
    auto &producer = wells[1];
    producer.control = MPMC::WellControl::ReservoirTotalRate;
    producer.primaryControl = MPMC::WellControl::ReservoirTotalRate;
    return wells;
}

void initializeSolution(Grid &grid, Vec solution)
{
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(solution, 0.0));
    auto values =
        grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();

    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                auto &primary = values[k][j][i];
                primary.fill(0.0);
                primary[Indices::Primary::pressure] =
                    CaseConfig::InitialState::pressure;
                for (int component = 0;
                     component < Indices::numIndependentCompositionsPerPhase;
                     ++component)
                {
                    const auto c = static_cast<std::size_t>(component);
                    primary[Indices::Primary::liquidComposition[c]] =
                        CaseConfig::InitialState::oilComposition[c];
                    primary[Indices::Primary::vaporComposition[c]] =
                        CaseConfig::InitialState::gasComposition[c];
                }
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

    const double scwFraction = readScwInjectionFraction();
    auto run = MPMC::cases::readRunOptions<Config>();
    char requested[PETSC_MAX_PATH_LEN]{};
    PetscBool resultSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetString(
        nullptr, nullptr, "-result_dir", requested, sizeof(requested),
        &resultSet));
    if (!resultSet)
        run.resultDirectory = scwFraction <= 1.0e-12
            ? "./results/co2-only"
            : "./results/scw-co2";
    MPMC::cases::printRunSummary(CaseConfig::name, run, rank);

    Grid grid(CaseConfig::Grid::nx, CaseConfig::Grid::ny,
              CaseConfig::Grid::nz,
              MPMC::GridExtent{CaseConfig::Grid::lx,
                               CaseConfig::Grid::ly,
                               CaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = makeFluid();
    const double waterDensity = MPMC::IapwsIf97WaterDensity::density(
        ScwMigration3D::initialPressure, ScwMigration3D::temperature);
    const double waterViscosity = MPMC::scalarValue(fluid.waterViscosity(
        ScwMigration3D::initialPressure));
    if (rank == 0)
        PetscPrintf(PETSC_COMM_SELF,
            "[BENCHMARK] model=independent-SCW+PR(CO2-nC4-nC16) "
            "grid=24x16x6 T=%.8g K P=%.8g Pa nominalPV=%.12g m3\n"
            "[SCW] IAPWS density=%.12g kg/m3 McBride-Wright viscosity=%.12g Pa.s\n"
            "[INJECTION] reservoir-total=%.12g m3/s SCW/CO2 volume fraction=%.8g/%.8g\n"
            "[INIT] S(O/G/W)=%.8g/%.8g/%.8g oil-z(CO2/nC4/nC16)=%.8g/%.8g/%.8g\n",
            ScwMigration3D::temperature, ScwMigration3D::initialPressure,
            ScwMigration3D::nominalPoreVolume, waterDensity,
            waterViscosity, ScwMigration3D::totalInjectionRate,
            scwFraction, 1.0 - scwFraction,
            CaseConfig::InitialState::oilSaturation,
            CaseConfig::InitialState::gasSaturation,
            CaseConfig::InitialState::waterSaturation,
            CaseConfig::InitialState::oilComposition[0],
            CaseConfig::InitialState::oilComposition[1],
            CaseConfig::InitialState::oilComposition[2]);

    auto runtimeOptions =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);
    runtimeOptions.scaling.enabled = true;
    runtimeOptions.scaling.rateWellResidualFloor =
        ScwMigration3D::totalInjectionRate;
    runtimeOptions.scaling.validate();

    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid, scwFraction));
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
