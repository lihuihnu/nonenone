/**
 * @file 3p6c_original.cpp
 * @brief 三相六组分原始对照算例的可执行程序入口与初始化流程。
 */
#include <case/petsc_custom_hooks.hpp>
#include <case/petsc_case_main.hpp>
#include <case/legacy_runtime_initialization.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_support.hpp>

#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Case
{

using Config = CaseConfig::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;



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
                    CaseConfig::Rock::kx,
                    CaseConfig::Rock::ky,
                    CaseConfig::Rock::kz};
                rock.porosity(i, j, k) = CaseConfig::Rock::porosity;
            }
}

double wellIndex(Grid &grid, int i, int j, int k,
                 double radius, double skin)
{
    const auto size = grid.cellSize({i, j, k});
    return MPMC::verticalPeacemanWellIndex(
        {size[0], size[1], size[2],
         CaseConfig::Rock::kx, CaseConfig::Rock::ky,
         radius, skin});
}

std::vector<Well> makeWells(Grid &grid)
{
    return MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid,
        WellConfig::wells,
        [](Grid &wellGrid, const auto &definition, int i, int j, int k)
        {
            return wellIndex(
                wellGrid, i, j, k, definition.radius, definition.skin);
        });
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
                auto &x = values[k][j][i];
                x.fill(0.0);
                x[Indices::Primary::pressure] = CaseConfig::InitialState::pressure;
                for (int c = 0; c < Indices::numIndependentCompositionsPerPhase; ++c)
                {
                    const auto n = static_cast<std::size_t>(c);
                    x[Indices::Primary::liquidComposition[n]] = CaseConfig::InitialState::oilComposition[n];
                    x[Indices::Primary::vaporComposition[n]] = CaseConfig::InitialState::gasComposition[n];
                }
                if constexpr (Indices::hasWater)
                    x[Indices::Primary::waterSaturation] = CaseConfig::InitialState::waterSaturation;
                x[Indices::Primary::liquidSaturation] = CaseConfig::InitialState::oilSaturation;
                x[Indices::Primary::vaporSaturation] = CaseConfig::InitialState::gasSaturation;

                if constexpr (Indices::hasAqueousCO2Dissolution)
                    x[Indices::Primary::aqueousCO2MoleFraction] =
                        CaseConfig::Dissolution::initialWaterCO2MoleFraction;
            }

    // 井 BHP 初值只写代表单元。
    if constexpr (Indices::hasWellUnknown)
    {
        for (const auto &def : WellConfig::wells)
        {
            const int i = MPMC::cases::resolveStructuredIndex(def.completion.i, grid.dimensions()[0]);
            const int j = MPMC::cases::resolveStructuredIndex(def.completion.j, grid.dimensions()[1]);
            const int k = def.completion.kBegin;
            if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
                j >= owned.yStart && j < owned.yStart + owned.yCount &&
                k >= owned.zStart && k < owned.zStart + owned.zCount)
                values[k][j][i][Indices::Primary::wellPressure] = def.initialBhp;
        }
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPMC::cases::validateCaseConfig<Indices, Config>();
    const auto run = MPMC::cases::readRunOptions<Config>();
    MPMC::cases::printRunSummary(CaseConfig::name, run, rank);

    Grid grid(CaseConfig::Grid::nx, CaseConfig::Grid::ny, CaseConfig::Grid::nz,
              MPMC::GridExtent{CaseConfig::Grid::lx, CaseConfig::Grid::ly, CaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = MPMC::cases::makeFluidSystem<Indices, Config>();
    Runtime runtime(grid, fluid,
                    MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run));
    runtime.setWells(makeWells(grid));

    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, solution);
    MPMC::cases::initializeLegacyRuntime<
        Indices, CaseConfig::InitialState>(runtime, solution);

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
    return MPMC::cases::runPetscCaseMain(argc, argv, [] { return Case::run(); });
}
