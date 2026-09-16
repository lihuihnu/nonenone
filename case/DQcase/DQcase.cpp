/**
 * @file DQcase.cpp
 * @brief DQ 网格三相多组分超算算例的可执行程序入口与初始化流程。
 */
#include <case/petsc_custom_hooks.hpp>
#include <case/petsc_case_main.hpp>
#include <case/legacy_runtime_initialization.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_support.hpp>
#include <common/console.hpp>

#include <cpgrid/cpgrid.hpp>
#include <cpgrid/grdecl.hpp>
#include <cpgrid/distributed_mesh_loader.hpp>
#include <cpgrid/grid_report.hpp>
#include <cpgrid/mesh.hpp>
#include <cpgrid/rock_loader.hpp>
#include <cpgrid/distributed_rock_loader.hpp>
#include <indices/indices.hpp>
#include <natural/petsc/natural_cpgrid_runtime.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
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
using Grid = MPMC::CpGridCore;
using Runtime = MPMC::NaturalCpGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;


std::vector<Well> makeWells(MPMC::Mesh &mesh)
{
    if constexpr (!Indices::hasWellUnknown)
        return {};

    std::vector<Well> result;
    result.reserve(WellConfig::wells.size());

    for (const auto &definition : WellConfig::wells)
    {
        if (definition.perforations.empty())
            throw std::invalid_argument(
                std::string(definition.name) + " has no perforation.");

        std::vector<MPMC::WellPerforation<PetscInt>> perforations;
        perforations.reserve(definition.perforations.size());
        for (const auto &perforation : definition.perforations)
        {
            const auto currentCell = static_cast<PetscInt>(
                mesh.currentIdFromInputIndex(perforation.inputCell));
            perforations.push_back({currentCell, perforation.wellIndex});
        }

        const PetscInt representativeCell = perforations.front().currentCellId;
        result.push_back(MPMC::cases::makeWell<Indices, PetscInt>(
            definition, representativeCell, std::move(perforations)));
    }
    return result;
}

void initializeSolution(Grid &grid, MPMC::Mesh &mesh, Vec solution)
{
    PetscCallAbort(mesh.communicator(), VecSet(solution, 0.0));
    const auto &map = grid.dofMap(Indices::numPrimaryVariables);

    PetscInt begin = 0, end = 0;
    PetscScalar *array = nullptr;
    PetscCallAbort(mesh.communicator(), VecGetOwnershipRange(solution, &begin, &end));
    PetscCallAbort(mesh.communicator(), VecGetArray(solution, &array));

    std::unordered_map<PetscInt, double> initialBhp;
    if constexpr (Indices::hasWellUnknown)
    {
        for (const auto &def : WellConfig::wells)
        {
            if (!def.perforations.empty())
                initialBhp.emplace(
                    mesh.currentIdFromInputIndex(def.perforations.front().inputCell),
                    def.initialBhp);
        }
    }

    for (PetscInt currentId : mesh.ownedCellIds())
    {
        const auto &cell = mesh.cellByCurrentId(currentId);
        std::array<double, Indices::numPrimaryVariables> x{};
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

        if constexpr (Indices::hasWellUnknown)
        {
            const auto it = initialBhp.find(cell.id());
            if (it != initialBhp.end())
                x[Indices::Primary::wellPressure] = it->second;
        }

        for (int v = 0; v < Indices::numPrimaryVariables; ++v)
        {
            const PetscInt global = map.globalIndex(cell, v);
            if (global < begin || global >= end)
                throw std::runtime_error("CpGrid owned-cell DOF is outside Vec ownership range.");
            array[global - begin] = x[static_cast<std::size_t>(v)];
        }
    }

    PetscCallAbort(mesh.communicator(), VecRestoreArray(solution, &array));
}

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPMC::cases::validateCaseConfig<Indices, Config>();
    auto run = MPMC::cases::readRunOptions<Config>();

    char grdeclPath[PETSC_MAX_PATH_LEN]{};
    PetscBool grdeclSet = PETSC_FALSE;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetString(
            nullptr, nullptr, "-grdecl",
            grdeclPath, sizeof(grdeclPath), &grdeclSet));

    if (grdeclSet == PETSC_TRUE)
        run.meshDirectory = grdeclPath;
    if (run.meshDirectory.empty())
        throw std::invalid_argument(
            "DQcase requires -grdecl, -mesh_dir, or CaseConfig::Grid::meshDirectory.");
    MPMC::cases::printRunSummary(CaseConfig::name, run, rank);

    MPMC::GrdeclLoadOptions loadOptions;
    // 与生成 DQ_data 的 MRST buildGridFromGRDECL(..., "largest") 保持一致，
    // 避免孤立 ACTNUM 岛形成与主储层无关的奇异线性子系统。
    loadOptions.keepLargestConnectedComponent = true;

    auto distributedInput = MPMC::loadDistributedMeshFromRoot(
        run.meshDirectory,
        grdeclSet == PETSC_TRUE,
        loadOptions,
        PETSC_COMM_WORLD,
        0);
    MPMC::Mesh &mesh = *distributedInput.mesh;

    // External numbering is always the original zero-based grid-file row.
    // Save the permutation explicitly so a result row can be audited against the
    // partition/current-id ordering used internally by PETSc.
    if (rank == 0)
        std::filesystem::create_directories(run.resultDirectory);
    MPI_Barrier(PETSC_COMM_WORLD);
    mesh.writeCellIdMap(run.resultDirectory + "/cell_id_map.csv");
    if (rank == 0)
    {
        MPMC::ConsoleSection indexing("GRID NUMBERING");
        indexing.row("External cell id", "input/active-cell row")
            .row("PETSc current id", "internal partition ordering")
            .row("Cell id map", run.resultDirectory + "/cell_id_map.csv");
        const std::string text = indexing.str();
        PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
    }

    if (rank == 0 && distributedInput.grdeclOnRoot &&
        distributedInput.grdeclOnRoot->disconnectedCellCountRemoved > 0)
    {
        MPMC::ConsoleSection connectivity("GRDECL CONNECTIVITY");
        connectivity.row(
            "Disconnected cells removed",
            std::to_string(distributedInput.grdeclOnRoot->disconnectedCellCountRemoved))
            .row("Component policy", "largest Cartesian face-connected component");
        const std::string text = connectivity.str();
        PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
    }

    std::optional<MPMC::CpGridRootRockData> rootGrdeclRock;
    if (distributedInput.format == MPMC::DistributedMeshInputFormat::Grdecl &&
        rank == 0 && distributedInput.grdeclOnRoot)
    {
        rootGrdeclRock = MPMC::takeRootRockData(*distributedInput.grdeclOnRoot);
        // Mesh 已独立拥有拓扑；只保留 PORO/PERM，尽早释放 expanded corners。
        distributedInput.grdeclOnRoot.reset();
    }

    /*
     * 全局编号表已经写出，后续求解只需要 owned+one-ring ghost 拓扑。
     * non-root 在 root-only ingest 阶段已经是 local snapshot；这里让 root
     * 同样释放远端 Node/Face/Polyhedron，避免 rank 0 长期承担 O(Nglobal)
     * 完整 Mesh 内存。
     */
    mesh.compactToLocalSnapshot(false);

    Grid grid(mesh);

    MPMC::CpGridRockLoadOptions rockOptions;
    rockOptions.replaceZeroPermeabilityWithMean =
        CaseConfig::Grid::replaceZeroPermeabilityWithMean;
    rockOptions.replaceZeroPorosityWithMean =
        CaseConfig::Grid::replaceZeroPorosityWithMean;
    if (distributedInput.format == MPMC::DistributedMeshInputFormat::Grdecl)
    {
        // 数值：生成参考 DQ_data 时，grid_created.m 会把零岩石属性替换为
        // 正值样本均值；仅在 direct GRDECL 路径复现该预处理，CSV 已完成归一化。
        rockOptions.replaceZeroPermeabilityWithMean = true;
        rockOptions.replaceZeroPorosityWithMean = true;
        MPMC::loadRockPropertiesFromRootData(
            grid,
            rootGrdeclRock ? &*rootGrdeclRock : nullptr,
            rockOptions,
            0);
    }
    else
    {
        MPMC::loadRockPropertiesFromCsvRoot(
            grid,
            run.meshDirectory,
            rockOptions,
            0);
    }

    // PETSc Vec 已完成装配，root 不再保留完整 PORO/PERM 输入数组。
    rootGrdeclRock.reset();
    // CSV 路径没有 grdeclOnRoot；GRDECL 路径已在 rock 装配前释放其 geometry。
    distributedInput.grdeclOnRoot.reset();
    grid.setup();
    MPMC::printCpGridSummary(grid);

    auto fluid = MPMC::cases::makeFluidSystem<Indices, Config>();
    Runtime runtime(grid, fluid,
                    MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run));
    runtime.setWells(makeWells(mesh));

    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, mesh, solution);
    MPMC::cases::initializeLegacyRuntime<
        Indices, CaseConfig::InitialState>(runtime, solution);

    MPMC::cases::NaturalSolver<Runtime> solver(runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, Config>(
        runtime, solver.snes(), solution, run, rank);

    PetscCallAbort(mesh.communicator(), VecDestroy(&solution));
    return 0;
}

} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return MPMC::cases::runPetscCaseMain(argc, argv, [] { return Case::run(); });
}
