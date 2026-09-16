/**
 * @file structured_three_eos_compare_runner.hpp
 * @brief StructuredGrid 全组分 PR/SW/CPA 对比算例的公共运行器。
 *
 * Case-specific physics stays in each case_config.hpp / well_config.hpp.  A
 * Definition only supplies the common Config type, the three EOS factory
 * configs, and the case's well table.  The runner owns only the duplicated
 * executable orchestration shared by the comparison cases.
 */
#pragma once

#include <case/case_support.hpp>
#include <case/natural_scaling_options.hpp>
#include <case/petsc_case_main.hpp>
#include <case/well_factory.hpp>

#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MPMC::cases
{

/**
 * @brief 统一 StructuredGrid PR/SW/CPA 对比算例的可执行运行流程。
 *
 * Required Definition contract:
 * - `using Config = ...;`
 * - `using PrFactoryConfig = ...;`
 * - `using SwFactoryConfig = ...;`
 * - `using CpaFactoryConfig = ...;`
 * - `static constexpr const auto &wells() noexcept;`
 *
 * The Config follows the normal case_config.hpp contract.  No grid topology,
 * rock parameter, well parameter, EOS parameter, solver setting, or output
 * policy is stored here.
 */
template <class Definition>
class StructuredThreeEosCompareRunner final
{
public:
    using Config = typename Definition::Config;
    using ModelConfig = MPMC::CompositionalModelConfig<
        Config::Model::numberOfComponents,
        Config::Model::hasWater,
        Config::Model::hasWells,
        Config::Model::enableDissolution,
        Config::Model::enableAdsorption,
        Config::Model::enableLandTrapping,
        Config::Model::phaseBehavior>;
    using Indices = MPMC::ADIndices<ModelConfig>;
    using Grid = MPMC::StructuredGridCore;
    using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
    using Well = MPMC::NaturalWell<Indices>;

    static_assert(Indices::fullyCompositionalThreePhase,
                  "This case must use the fully compositional O/G/W model.");

    enum class EosChoice { Pr, Sw, Cpa };

    [[nodiscard]] static const char *eosToken(EosChoice choice)
    {
        switch (choice)
        {
        case EosChoice::Pr:  return "pr";
        case EosChoice::Sw:  return "sw";
        case EosChoice::Cpa: return "cpa";
        }
        return "unknown";
    }

    [[nodiscard]] static const char *eosDisplayName(EosChoice choice)
    {
        switch (choice)
        {
        case EosChoice::Pr:  return "Peng-Robinson";
        case EosChoice::Sw:  return "Soreide-Whitson";
        case EosChoice::Cpa: return "Cubic-Plus-Association";
        }
        return "unknown";
    }

    [[nodiscard]] static EosChoice readEosChoice()
    {
        char value[32] = "pr";
        PetscBool set = PETSC_FALSE;
        PetscCallAbort(PETSC_COMM_WORLD,
                       PetscOptionsGetString(nullptr, nullptr, "-eos",
                                             value, sizeof(value), &set));
        const std::string token(value);
        if (token == "pr" || token == "PR")
            return EosChoice::Pr;
        if (token == "sw" || token == "SW" || token == "soreide-whitson")
            return EosChoice::Sw;
        if (token == "cpa" || token == "CPA")
            return EosChoice::Cpa;
        throw std::invalid_argument("-eos must be one of: pr, sw, cpa.");
    }

    [[nodiscard]] static MPMC::FluidSystem<Indices> makeComparisonFluid(EosChoice choice)
    {
        switch (choice)
        {
        case EosChoice::Pr:
            return MPMC::cases::makeFluidSystem<Indices, typename Definition::PrFactoryConfig>();
        case EosChoice::Sw:
            return MPMC::cases::makeFluidSystem<Indices, typename Definition::SwFactoryConfig>();
        case EosChoice::Cpa:
            return MPMC::cases::makeFluidSystem<Indices, typename Definition::CpaFactoryConfig>();
        }
        throw std::logic_error("Unknown EOS choice.");
    }

    static void printInitialFlash(
        const MPMC::FluidSystem<Indices> &fluid,
        EosChoice choice,
        PetscMPIInt rank)
    {
        if (rank != 0)
            return;
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = Config::Fluid::waterComponent;
        const MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
        const auto result = flash.flash(
            Config::InitialState::pressure,
            Config::InitialState::temperature,
            Config::InitialState::overallComposition);
        if (!result.converged)
            throw std::runtime_error("Selected EOS cannot flash the common initial P-T-z state.");

        PetscPrintf(PETSC_COMM_SELF,
                    "[EOS] %s (%s)\n"
                    "[EOS] initial beta(O/G/W)=%.8g / %.8g / %.8g\n"
                    "[EOS] initial S(O/G/W)   =%.8g / %.8g / %.8g\n",
                    eosDisplayName(choice), eosToken(choice),
                    result.phaseMoleFraction[0], result.phaseMoleFraction[1], result.phaseMoleFraction[2],
                    result.saturation[0], result.saturation[1], result.saturation[2]);
    }

    static void initializeRock(Grid &grid)
    {
        grid.initializeRockProperties();
        auto rock = grid.rockWriteView();
        const auto owned = grid.ownedRegion();

        for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
            for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
                for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
                {
                    rock.permeability(i, j, k) = {
                        Config::Rock::kx(i, j, k),
                        Config::Rock::ky(i, j, k),
                        Config::Rock::kz(i, j, k)};
                    rock.porosity(i, j, k) = Config::Rock::porosity(i, j, k);
                }
    }

    [[nodiscard]] static double wellIndex(
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
             Config::Rock::kx(i, j, k), Config::Rock::ky(i, j, k),
             radius, skin});
    }

    [[nodiscard]] static std::vector<Well> makeWells(Grid &grid)
    {
        return MPMC::cases::makeStructuredWells<Indices, PetscInt>(
            grid,
            Definition::wells(),
            [](Grid &wellGrid, const auto &definition, int i, int j, int k)
            {
                return wellIndex(
                    wellGrid, i, j, k, definition.radius, definition.skin);
            });
    }

    /**
     * @brief 在统一 P-T-z flash 初始化后覆盖代表单元的 BHP 初值。
     *
     * `initializeUniformFromPTZ()` initializes all reservoir variables and first
     * sets well-pressure DOFs to reservoir pressure.  This changes only the
     * representative well-cell BHP initial guess, not the fluid state.
     */
    static void initializeWellPressureGuesses(Grid &grid, Vec solution)
    {
        auto values = grid.template vecGetArray<Indices::numPrimaryVariables>(solution);
        const auto owned = grid.ownedRegion();
        const auto dims = grid.dimensions();

        for (const auto &def : Definition::wells())
        {
            const int i = MPMC::cases::resolveStructuredIndex(def.completion.i, dims[0]);
            const int j = MPMC::cases::resolveStructuredIndex(def.completion.j, dims[1]);
            const int k = def.completion.kBegin;
            if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
                j >= owned.yStart && j < owned.yStart + owned.yCount &&
                k >= owned.zStart && k < owned.zStart + owned.zCount)
            {
                values[k][j][i][Indices::Primary::wellPressure] = def.initialBhp;
            }
        }
        grid.template vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
    }

    static int run()
    {
        PetscMPIInt rank = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        MPMC::cases::validateCaseConfig<Indices, Config>();
        const EosChoice eosChoice = readEosChoice();
        auto run = MPMC::cases::readRunOptions<Config>();

        // Default to one result directory per EOS; explicit -result_dir wins.
        char requestedResult[PETSC_MAX_PATH_LEN]{};
        PetscBool resultSet = PETSC_FALSE;
        PetscCallAbort(PETSC_COMM_WORLD,
                       PetscOptionsGetString(nullptr, nullptr, "-result_dir",
                                             requestedResult, sizeof(requestedResult), &resultSet));
        if (resultSet == PETSC_FALSE)
            run.resultDirectory = std::string("./results/") + eosToken(eosChoice);

        MPMC::cases::printRunSummary(Config::name, run, rank);

        Grid grid(
            Config::Grid::nx,
            Config::Grid::ny,
            Config::Grid::nz,
            MPMC::GridExtent{
                Config::Grid::lx,
                Config::Grid::ly,
                Config::Grid::lz});
        grid.setup();
        initializeRock(grid);
        MPMC::printStructuredGridSummary(grid);

        auto fluid = makeComparisonFluid(eosChoice);
        printInitialFlash(fluid, eosChoice, rank);

        auto runtimeOptions =
            MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);

        // Natural scaling 只改变方程坐标，不改变物理零点。
        MPMC::cases::applyNaturalScalingCaseDefaults<typename Config::Numerics>(
            runtimeOptions);
        MPMC::cases::applyNaturalScalingPetscOptions(runtimeOptions, rank);

        Runtime runtime(grid, fluid, runtimeOptions);
        runtime.setWells(makeWells(grid));

        Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);

        // Fully-compositional initialization is defined only by P, T, and z;
        // flash supplies phase compositions, fractions, Z factors, and saturation.
        runtime.initializeUniformFromPTZ(
            solution,
            Config::InitialState::pressure,
            Config::InitialState::temperature,
            Config::InitialState::overallComposition);
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

    static int runPetscMain(int argc, char **argv)
    {
        return MPMC::cases::runPetscCaseMain(argc, argv, [] { return run(); });
    }
};

} // namespace MPMC::cases
