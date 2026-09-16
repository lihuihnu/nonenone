/**
 * @file panfili2025_case3_fullphysics.cpp
 * @brief Panfili 2025 Case-3 全物理算例的可执行程序入口与初始化流程。
 */
#include <case/petsc_custom_hooks.hpp>
#include <case/petsc_case_main.hpp>
#include "case_config.hpp"
#include "case_fluid.hpp"
#include "well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_support.hpp>

#include <common/math.hpp>
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
    CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;
using Eval = typename Indices::ValueType;

static_assert(Indices::fullyCompositionalThreePhase,
              "Panfili-2025 case #3 requires the live-water O/G/W formulation.");



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

double wellIndex(Grid &grid, int i, int j, int k, double radius, double skin)
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

void setPrimaryOperatingMode(
    Well &well,
    MPMC::WellType type,
    MPMC::WellControl control,
    double target)
{
    well.type = type;
    well.control = control;
    well.target = MPMC::isRateControl(control) ? std::abs(target) : target;
    well.primaryControl = control;
    well.primaryTarget = MPMC::isRateControl(control) ? std::abs(target) : target;
    well.primaryControlInitialized = true;
}

Runtime::WellScheduleUpdater makePanfiliWellScheduleUpdater()
{
    // State is kept inside the updater so ordinary BHP/rate limit switching is
    // NOT reset on every time step.  Primary control is reset only when the
    // physical operating stage changes (depletion / idle / injection / closed).
    return [stage = std::array<WellConfig::OperatingStage, 2>{
                WellConfig::OperatingStage::Closed,
                WellConfig::OperatingStage::Closed}](
               double time,
               std::vector<Well> &wells) mutable
    {
        if (wells.size() != WellConfig::wells.size())
            throw std::logic_error(
                "Panfili well schedule expects exactly two physical crest wells.");

        const auto desired = WellConfig::stageAt(time);
        for (std::size_t i = 0; i < wells.size(); ++i)
        {
            auto &well = wells[i];
            if (stage[i] == desired)
                continue;

            stage[i] = desired;
            switch (desired)
            {
            case WellConfig::OperatingStage::Depletion:
                well.schedule.enabled = true;
                setPrimaryOperatingMode(
                    well,
                    MPMC::WellType::Producer,
                    MPMC::WellControl::GasRate,
                    WellConfig::producerTarget());
                break;

            case WellConfig::OperatingStage::Idle:
                well.schedule.enabled = false;
                break;

            case WellConfig::OperatingStage::Injection:
                well.schedule.enabled = true;
                setPrimaryOperatingMode(
                    well,
                    MPMC::WellType::Injector,
                    MPMC::WellControl::GasRate,
                    WellConfig::injectorTarget());
                break;

            case WellConfig::OperatingStage::Closed:
                well.schedule.enabled = false;
                break;
            }
        }
    };
}

void initializeWellPressureGuesses(Grid &grid, Vec solution)
{
    auto values = grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();
    const auto dims = grid.dimensions();

    for (const auto &def : WellConfig::wells)
    {
        const int i = MPMC::cases::resolveStructuredIndex(def.completion.i, dims[0]);
        const int j = MPMC::cases::resolveStructuredIndex(def.completion.j, dims[1]);
        const int k = def.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
        {
            values[k][j][i][Indices::Primary::wellPressure] =
                CaseConfig::InitialState::pressure;
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

    Grid grid(
        CaseConfig::Grid::nx,
        CaseConfig::Grid::ny,
        CaseConfig::Grid::nz,
        MPMC::GridExtent{
            CaseConfig::Grid::lx,
            CaseConfig::Grid::ly,
            CaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = MPMC::cases::makeFluidSystem<Indices, Config>();
    Panfili2025CaseFluid::applyRockFluidApproximation(fluid);

    Runtime runtime(
        grid,
        fluid,
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run));
    runtime.setWells(makeWells(grid));
    runtime.setWellScheduleUpdater(makePanfiliWellScheduleUpdater());

    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    runtime.initializeUniformFromPTZ(
        solution,
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
    initializeWellPressureGuesses(grid, solution);
    runtime.initializeHistory(solution);

    MPMC::cases::NaturalSolver<Runtime> solver(
        runtime,
        grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, Config>(
        runtime,
        solver.snes(),
        solution,
        run,
        rank);

    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));
    return 0;
}

} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return MPMC::cases::runPetscCaseMain(argc, argv, [] { return Case::run(); });
}
