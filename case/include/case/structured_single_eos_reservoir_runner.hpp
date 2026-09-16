/**
 * @file structured_single_eos_reservoir_runner.hpp
 * @brief StructuredGrid 单 EOS 全组分三相储层算例的公共运行器。
 *
 * 各算例继续在自己的 case_config.hpp / well_config.hpp 中保存全部物理参数。
 * Definition 只向公共运行流程暴露 Config 与井表；本文件不保存 EOS、网格、
 * 岩石、井、求解器或输出参数。
 */
#pragma once

#include <case/case_support.hpp>
#include <case/natural_scaling_options.hpp>
#include <case/petsc_case_main.hpp>
#include <case/well_factory.hpp>

#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <type_traits>
#include <vector>

namespace MPMC::cases
{

namespace detail
{

template <class Numerics, class = void>
struct HasNaturalScalingConfiguration : std::false_type {};

template <class Numerics>
struct HasNaturalScalingConfiguration<Numerics, std::void_t<
    decltype(Numerics::pressureScale),
    decltype(Numerics::compositionScale),
    decltype(Numerics::saturationScale),
    decltype(Numerics::massResidualScale),
    decltype(Numerics::fugacityResidualScale),
    decltype(Numerics::closureResidualScale),
    decltype(Numerics::rateWellResidualFloor)>> : std::true_type {};

} // namespace detail

/**
 * @brief 统一单 EOS StructuredGrid 全组分三相算例的可执行运行流程。
 *
 * Definition 需要提供：
 * - `using Config = ...;`
 * - `static constexpr const auto &wells() noexcept;`
 *
 * Config 继续遵循现有 case_config.hpp 契约。公共运行器只承载两个算例原先
 * 完全重复的启动、初始化、求解与收尾流程，不改变任何物理或数值配置。
 */
template <class Definition>
class StructuredSingleEosReservoirRunner final
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
                        Config::Rock::kx,
                        Config::Rock::ky,
                        Config::Rock::kz};
                    rock.porosity(i, j, k) = Config::Rock::porosity;
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
             Config::Rock::kx, Config::Rock::ky,
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
     * `initializeUniformFromPTZ()` 先从热力学初始化全部储层变量，并把井压
     * 自由度设为储层压力作为中性初值；这里只覆盖代表井单元的 BHP 初猜，
     * 不改变流体状态输入。
     */
    static void initializeWellPressureGuesses(Grid &grid, Vec solution)
    {
        auto values =
            grid.template vecGetArray<Indices::numPrimaryVariables>(solution);
        const auto owned = grid.ownedRegion();
        const auto dims = grid.dimensions();

        for (const auto &def : Definition::wells())
        {
            const int i = MPMC::cases::resolveStructuredIndex(
                def.completion.i, dims[0]);
            const int j = MPMC::cases::resolveStructuredIndex(
                def.completion.j, dims[1]);
            const int k = def.completion.kBegin;
            if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
                j >= owned.yStart && j < owned.yStart + owned.yCount &&
                k >= owned.zStart && k < owned.zStart + owned.zCount)
            {
                values[k][j][i][Indices::Primary::wellPressure] =
                    def.initialBhp;
            }
        }
        grid.template vecRestoreArray<Indices::numPrimaryVariables>(
            solution, values);
    }

    static int run()
    {
        PetscMPIInt rank = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        MPMC::cases::validateCaseConfig<Indices, Config>();
        const auto run = MPMC::cases::readRunOptions<Config>();
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

        auto fluid = MPMC::cases::makeFluidSystem<Indices, Config>();
        auto runtimeOptions =
            MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);
        if constexpr (detail::HasNaturalScalingConfiguration<
                          typename Config::Numerics>::value)
        {
            runtimeOptions.scaling.pressureScale =
                Config::Numerics::pressureScale;
            runtimeOptions.scaling.compositionScale =
                Config::Numerics::compositionScale;
            runtimeOptions.scaling.saturationScale =
                Config::Numerics::saturationScale;
            runtimeOptions.scaling.massResidualScale =
                Config::Numerics::massResidualScale;
            runtimeOptions.scaling.fugacityResidualScale =
                Config::Numerics::fugacityResidualScale;
            runtimeOptions.scaling.closureResidualScale =
                Config::Numerics::closureResidualScale;
            runtimeOptions.scaling.rateWellResidualFloor =
                Config::Numerics::rateWellResidualFloor;
            MPMC::cases::applyNaturalScalingPetscOptions(runtimeOptions, rank);
        }
        Runtime runtime(grid, fluid, runtimeOptions);
        runtime.setWells(makeWells(grid));

        Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);

        // 全组分初始化只由 P、T 和总体组成 z 定义；相组成、相分率、Z 因子与
        // 饱和度均由 flash 产生，不作为额外用户输入。
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
