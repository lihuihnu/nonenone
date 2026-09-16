/**
 * @file five_component_eos_tuned_legacy_compare.cpp
 * @brief PR 干油气 Flash 与独立不可混相水相耦合的传统对照算例。
 */
#include <case/petsc_custom_hooks.hpp>
#include <case/petsc_case_main.hpp>
#include <case/natural_scaling_options.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_support.hpp>

#include <indices/indices.hpp>
#include <natural/phase_state.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace Case
{

using Config = LegacyCaseConfig::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    LegacyCaseConfig::Model::numberOfComponents,
    LegacyCaseConfig::Model::hasWater,
    LegacyCaseConfig::Model::hasWells,
    LegacyCaseConfig::Model::enableDissolution,
    LegacyCaseConfig::Model::enableAdsorption,
    LegacyCaseConfig::Model::enableLandTrapping,
    LegacyCaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;

static_assert(!Indices::fullyCompositionalThreePhase,
              "This case must use the legacy independent-water model.");
static_assert(Indices::hasIndependentWaterConservation,
              "This case requires an independent water conservation equation.");

MPMC::FluidSystem<Indices> makeComparisonFluid()
{
    return MPMC::cases::makeFluidSystem<
        Indices, LegacyCaseConfig::PrFactoryConfig>();
}

Flash::Result initialDryFlash(
    const MPMC::FluidSystem<Indices> &fluid,
    PetscMPIInt rank)
{
    MPMC::PhasePresence oilGas = MPMC::PhasePresence::oilOnly();
    oilGas.add(MPMC::CompositionalPhase::Gas);
    const Flash flash(fluid.eos);
    const auto result = flash.flashRestricted(
        LegacyCaseConfig::InitialState::pressure,
        LegacyCaseConfig::InitialState::temperature,
        LegacyCaseConfig::InitialState::dryOverallComposition,
        oilGas);
    if (!result.converged ||
        !result.presence.contains(MPMC::CompositionalPhase::Oil) ||
        !result.presence.contains(MPMC::CompositionalPhase::Gas))
    {
        throw std::runtime_error(
            "Selected EOS cannot produce the required dry oil+gas initial equilibrium.");
    }

    if (rank == 0)
    {
        PetscPrintf(PETSC_COMM_SELF,
                    "[EOS] %s (%s), legacy independent H2O\n"
                    "[EOS] dry beta(O/G)=%.8g / %.8g\n",
                    "Peng-Robinson dry oil/gas", "pr",
                    result.phaseMoleFraction[0], result.phaseMoleFraction[1]);
    }
    return result;
}

double matchedWaterSaturation(const Flash::Result &flash)
{
    const double dryMolarDensity =
        flash.saturation[0] * flash.molarDensity[0] +
        flash.saturation[1] * flash.molarDensity[1];
    const double waterMolarDensity =
        LegacyCaseConfig::CommonFluid::surfaceDensity[2] /
        LegacyCaseConfig::InitialState::waterMolarMass;
    const double zw =
        LegacyCaseConfig::InitialState::overallWaterMoleFraction;
    const double denominator =
        (1.0 - zw) * waterMolarDensity + zw * dryMolarDensity;
    const double sw = zw * dryMolarDensity / denominator;
    if (!(dryMolarDensity > 0.0) || !(waterMolarDensity > 0.0) ||
        !(sw > 0.0 && sw < 1.0) || !std::isfinite(sw))
        throw std::runtime_error(
            "Cannot reconstruct the requested overall H2O mole fraction.");
    return sw;
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
                    LegacyCaseConfig::Rock::kx(i, j, k),
                    LegacyCaseConfig::Rock::ky(i, j, k),
                    LegacyCaseConfig::Rock::kz(i, j, k)};
                rock.porosity(i, j, k) =
                    LegacyCaseConfig::Rock::porosity(i, j, k);
            }
}

double wellIndex(
    Grid &grid,
    int i,
    int j,
    int k,
    double radius,
    double skin)
{
    const auto size = grid.cellSize({i, j, k});
    return MPMC::verticalPeacemanWellIndex(
        {size[0], size[1], size[2],
         LegacyCaseConfig::Rock::kx(i, j, k),
         LegacyCaseConfig::Rock::ky(i, j, k),
         radius, skin});
}

std::vector<Well> makeWells(Grid &grid)
{
    return MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid,
        LegacyWellConfig::wells,
        [](Grid &wellGrid, const auto &definition, int i, int j, int k)
        {
            return wellIndex(
                wellGrid, i, j, k, definition.radius, definition.skin);
        });
}

void initializeSolution(
    Grid &grid,
    Vec solution,
    const Flash::Result &flash,
    double sw)
{
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(solution, 0.0));
    auto values = grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();
    const double hydrocarbon = 1.0 - sw;

    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                auto &primary = values[k][j][i];
                primary.fill(0.0);
                primary[Indices::Primary::pressure] =
                    LegacyCaseConfig::InitialState::pressure;
                for (int c = 0; c < Indices::numIndependentCompositionsPerPhase; ++c)
                {
                    const auto component = static_cast<std::size_t>(c);
                    primary[Indices::Primary::liquidComposition[component]] =
                        flash.composition[0][component];
                    primary[Indices::Primary::vaporComposition[component]] =
                        flash.composition[1][component];
                }
                primary[Indices::Primary::liquidSaturation] =
                    hydrocarbon * flash.saturation[0];
                primary[Indices::Primary::vaporSaturation] =
                    hydrocarbon * flash.saturation[1];
                primary[Indices::Primary::waterSaturation] = sw;
            }

    for (const auto &definition : LegacyWellConfig::wells)
    {
        const int i = MPMC::cases::resolveStructuredIndex(
            definition.completion.i, grid.dimensions()[0]);
        const int j = MPMC::cases::resolveStructuredIndex(
            definition.completion.j, grid.dimensions()[1]);
        const int k = definition.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
        {
            values[k][j][i][Indices::Primary::wellPressure] =
                definition.initialBhp;
        }
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    MPMC::cases::validateCaseConfig<Indices, Config>();
    auto run = MPMC::cases::readRunOptions<Config>();

    char requestedResult[PETSC_MAX_PATH_LEN]{};
    PetscBool resultSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetString(nullptr, nullptr, "-result_dir",
                                         requestedResult, sizeof(requestedResult), &resultSet));
    if (resultSet == PETSC_FALSE)
        run.resultDirectory = "./results/full/pr";

    MPMC::cases::printRunSummary(LegacyCaseConfig::name, run, rank);

    Grid grid(
        LegacyCaseConfig::Grid::nx,
        LegacyCaseConfig::Grid::ny,
        LegacyCaseConfig::Grid::nz,
        MPMC::GridExtent{
            LegacyCaseConfig::Grid::lx,
            LegacyCaseConfig::Grid::ly,
            LegacyCaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = makeComparisonFluid();
    const auto flash = initialDryFlash(fluid, rank);
    const double initialSw = matchedWaterSaturation(flash);
    if (rank == 0)
    {
        const double hydrocarbon = 1.0 - initialSw;
        PetscPrintf(PETSC_COMM_SELF,
                    "[EOS] initial S(O/G/W)=%.8g / %.8g / %.8g "
                    "(matched overall z_H2O=%.8g)\n",
                    hydrocarbon * flash.saturation[0],
                    hydrocarbon * flash.saturation[1], initialSw,
                    LegacyCaseConfig::InitialState::overallWaterMoleFraction);
    }
    auto runtimeOptions =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);

    MPMC::cases::applyNaturalScalingPetscOptions(runtimeOptions, rank);

    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid));

    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, solution, flash, initialSw);
    runtime.initializePhaseStateFromSolution(
        solution, MPMC::HydrocarbonPhaseState::TwoPhase);
    runtime.updateState(solution);
    runtime.initializeHistory(solution);

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
    return MPMC::cases::runPetscCaseMain(argc, argv, [] { return Case::run(); });
}
