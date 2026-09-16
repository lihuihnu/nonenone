/**
 * @file natural_cpgrid_test.cpp
 * @brief 集成测试：验证 `natural_cpgrid_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <cpgrid/cpgrid.hpp>
#include <cpgrid/distributed_mesh_loader.hpp>
#include <cpgrid/rock_loader.hpp>
#include <cpgrid/distributed_rock_loader.hpp>
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/petsc/natural_petsc.hpp>

#include <petscmat.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace AdapterTest
{

using ModelConfig =
    MPMC::CompositionalModelConfig<
        6,      // components
        true,   // water
        true,   // well unknown
        false,  // aqueous CO2 dissolution
        false,  // adsorption
        false>; // Land trapping

using Indices =
    MPMC::ADIndices<ModelConfig>;

using Grid =
    MPMC::CpGridCore;

using Runtime =
    MPMC::NaturalCpGridRuntime<
        Indices,
        Grid>;

struct DomainFailureRuntime final
{
    PetscErrorCode formFunction(SNES, Vec, Vec)
    {
        throw std::runtime_error("synthetic thermodynamic domain failure");
    }
};

template <class IndexType>
MPMC::CompositionalMixture<IndexType>
makeMixture()
{
    return MPMC::CompositionalMixture<IndexType>(
        {189.515,
         304.2,
         387.607,
         597.497,
         698.515,
         875.0},
        {4580011.59,
         7386592.50,
         4095515.97,
         3345244.875,
         1768374.5625,
         1169006.79},
        {9.97012032965401e-5,
         9.26344713533338e-5,
         0.000217076707259486,
         0.000381162235869935,
         0.000721410148917871,
         0.00113570073874421},
        {0.00854,
         0.228,
         0.16733,
         0.38609,
         0.80784,
         1.23141},
        {0.0161594,
         0.04401,
         0.0455725,
         0.11774,
         0.248827,
         0.48152},
        {{0.0,
          0.00070981,
          0.00077754,
          0.0100,
          0.0110,
          0.0110},
         {0.00070981,
          0.0,
          0.1500,
          0.1500,
          0.1500,
          0.1500},
         {0.00077754,
          0.1500,
          0.0,
          0.0,
          0.0,
          0.0},
         {0.0100,
          0.1500,
          0.0,
          0.0,
          0.0,
          0.0},
         {0.0110,
          0.1500,
          0.0,
          0.0,
          0.0,
          0.0},
         {0.0110,
          0.1500,
          0.0,
          0.0,
          0.0,
          0.0}});
}

MPMC::FluidSystem<Indices>
makeFluidSystem()
{
    MPMC::CubicEquationOfState<Indices> eos(
        0.4572355,
        0.0779691,
        makeMixture<Indices>(),
        1,
        2.414213562373095,
        -0.414213562373095);

    MPMC::FluidSystem<Indices> fluid(
        {800.0,
         2.0,
         1000.0},
        {1.0e-3,
         1.0e-5,
         2.0e-4},
        std::move(eos),
        {"N2/CH4",
         "CO2",
         "C2-5",
         "C6-13",
         "C14-24",
         "C25-80"},
        {"liquid",
         "vapor",
         "water"},
        387.45);

    using Eval = typename Indices::ValueType;

    fluid.gasRelativePermeability =
        [](Eval saturation)
        {
            return saturation * saturation;
        };

    fluid.waterRelativePermeability =
        [](Eval saturation)
        {
            return saturation * saturation;
        };

    fluid.threePhaseOilRelativePermeability =
        [](Eval,
           Eval oilSaturation,
           Eval)
        {
            return oilSaturation *
                   oilSaturation;
        };

    fluid.waterViscosity =
        [](Eval)
        {
            return Eval(2.0e-4);
        };

    fluid.waterFormationVolumeFactor =
        [](Eval)
        {
            return Eval(1.0);
        };

    return fluid;
}

void fillInitialSolution(
    Grid &grid,
    Vec solution)
{
    const auto &map =
        grid.dofMap(
            static_cast<PetscInt>(
                Indices::numPrimaryVariables));

    PetscInt begin = 0;
    PetscInt end = 0;
    PetscScalar *array = nullptr;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetOwnershipRange(
            solution,
            &begin,
            &end));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetArray(
            solution,
            &array));

    /*
     * 使用一组热力学有效、各单元略带压力梯度的状态。
     * 压力梯度避免有限差分时恰好落在 upstream tie-break 的非光滑点。
     */
    const std::array<double, 5>
        independentComposition{
            0.3246914,
            0.0128351,
            0.2278401,
            0.2606985,
            0.1134144};

    for (const auto &cell :
         grid.localCells())
    {
        std::array<double,
                   Indices::numPrimaryVariables>
            primary{};

        primary[
            Indices::Primary::pressure] =
            15.0e6 -
            100.0 *
                static_cast<double>(
                    cell.id());

        for (int component = 0;
             component <
                 Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            const std::size_t c =
                static_cast<std::size_t>(
                    component);

            primary[
                Indices::Primary::
                    liquidComposition[c]] =
                independentComposition[c];

            /*
             * 初值让两相组成相同，以避免测试本身依赖某个特定 flash 路径；
             * EOS 仍会分别取 liquid/vapor root。
             */
            primary[
                Indices::Primary::
                    vaporComposition[c]] =
                independentComposition[c];
        }

        primary[
            Indices::Primary::
                waterSaturation] = 0.20;

        primary[
            Indices::Primary::
                wellPressure] = 14.0e6;

        primary[
            Indices::Primary::
                liquidSaturation] = 0.50;

        primary[
            Indices::Primary::
                vaporSaturation] = 0.30;

        for (int component = 0;
             component <
                 Indices::numPrimaryVariables;
             ++component)
        {
            const PetscInt global =
                map.globalIndex(
                    cell,
                    component);

            if (global < begin ||
                global >= end)
            {
                throw std::runtime_error(
                    "Initial-state ownership is inconsistent with CpGrid DofMap.");
            }

            array[
                global - begin] =
                static_cast<PetscScalar>(
                    primary[
                        static_cast<std::size_t>(
                            component)]);
        }
    }

    PetscCallAbort(
        grid.mesh().communicator(),
        VecRestoreArray(
            solution,
            &array));
}

std::vector<MPMC::NaturalWell<Indices>>
makeCrossRankWell(const MPMC::Mesh &mesh)
{
    if (mesh.cellCount() < 2)
    {
        throw std::runtime_error(
            "Natural/CpGrid test requires at least two cells.");
    }

    MPMC::NaturalWell<Indices> well;
    well.type =
        MPMC::WellType::Producer;
    well.control =
        MPMC::WellControl::TotalRate;
    well.target = 0.0;
    well.bhpCellId = 0;

    well.perforations = {
        {0, 1.0e-12},
        {static_cast<PetscInt>(
             mesh.cellCount() - 1),
         1.0e-12}};

    return {well};
}

void requireFiniteVector(
    Vec vector,
    const char *label)
{
    PetscReal norm = 0.0;

    PetscCallAbort(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                vector)),
        VecNorm(
            vector,
            NORM_INFINITY,
            &norm));

    if (!std::isfinite(
            static_cast<double>(norm)))
    {
        throw std::runtime_error(
            std::string(label) +
            " contains a non-finite value.");
    }
}

/**
 * @brief 对一个 pressure 列做全局中心差分，验证完整 PETSc Jacobian。
 *
 * phase-state 在三次 FormFunction 之间保持冻结，这与单次 Newton 线性化一致。
 * 选 pressure 而不是 BHP，可避免 rate-control 的人为 -1 Jacobian shift
 * 对“严格 residual 导数”比较造成预期中的差异。
 */
double checkPressureJacobianColumn(
    SNES snes,
    Runtime &runtime,
    Vec solution,
    Mat jacobian)
{
    const auto &map =
        runtime.grid()
            .dofMap(
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables));

    const PetscInt pressureDof =
        map.globalIndex(
            PetscInt(0),
            Indices::Primary::pressure);

    Vec direction = nullptr;
    Vec plus = nullptr;
    Vec minus = nullptr;
    Vec fPlus = nullptr;
    Vec fMinus = nullptr;
    Vec finiteDifference = nullptr;
    Vec jacobianProduct = nullptr;

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &direction));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &plus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &minus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &fPlus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &fMinus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &finiteDifference));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDuplicate(
            solution,
            &jacobianProduct));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecSet(direction, 0.0));

    const PetscScalar one = 1.0;

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecSetValue(
            direction,
            pressureDof,
            one,
            INSERT_VALUES));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecAssemblyBegin(direction));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecAssemblyEnd(direction));

    constexpr double step = 1.0;

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecCopy(solution, plus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecCopy(solution, minus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecAXPY(plus, step, direction));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecAXPY(minus, -step, direction));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        SNESComputeFunction(
            snes,
            plus,
            fPlus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        SNESComputeFunction(
            snes,
            minus,
            fMinus));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecCopy(
            fPlus,
            finiteDifference));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecAXPY(
            finiteDifference,
            -1.0,
            fMinus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecScale(
            finiteDifference,
            0.5 / step));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        MatMult(
            jacobian,
            direction,
            jacobianProduct));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecAXPY(
            finiteDifference,
            -1.0,
            jacobianProduct));

    PetscReal differenceNorm = 0.0;
    PetscReal jacobianNorm = 0.0;

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecNorm(
            finiteDifference,
            NORM_INFINITY,
            &differenceNorm));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecNorm(
            jacobianProduct,
            NORM_INFINITY,
            &jacobianNorm));

    const double relativeError =
        static_cast<double>(differenceNorm) /
        std::max(
            1.0,
            static_cast<double>(jacobianNorm));

    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&direction));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&plus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&minus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&fPlus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&fMinus));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&finiteDifference));
    PetscCallAbort(
        runtime.grid().mesh().communicator(),
        VecDestroy(&jacobianProduct));

    return relativeError;
}

void runTest(const std::string &meshDirectory)
{
    auto distributed = MPMC::loadDistributedMeshFromRoot(
        meshDirectory,
        false,
        {},
        PETSC_COMM_WORLD,
        0);
    MPMC::Mesh &mesh = *distributed.mesh;

    // 与 DQ 生产路径一致：所有 rank 的求解拓扑都只保留 owned+ghost。
    mesh.compactToLocalSnapshot(false);

    Grid grid(mesh);

    MPMC::loadRockPropertiesFromCsvRoot(
        grid,
        meshDirectory,
        {},
        0);

    grid.setup();

    auto fluid = makeFluidSystem();

    Runtime::Options options;
    options.timeStep = 86400.0;
    options.fugacityScalingFactor = 1.0;

    Runtime runtime(
        grid,
        fluid,
        options);

    runtime.setWells(
        makeCrossRankWell(mesh));

    Vec solution =
        grid.createGlobalVector(
            static_cast<PetscInt>(
                Indices::numPrimaryVariables));

    fillInitialSolution(
        grid,
        solution);

    runtime.initializePhaseStateFromSolution(
        solution,
        MPMC::HydrocarbonPhaseState::TwoPhase);

    /*
     * 真实覆盖 updateState：归一化组成、phase appearance/disappearance、
     * K/z/L/Z secondary 更新。
     */
    runtime.updateState(solution);

    runtime.initializeHistory(solution);

    Vec residual =
        runtime.createResidualVector();

    Mat jacobian =
        runtime.createJacobian();

    // G8J: createJacobian() now materializes and final-assembles the complete
    // static block topology before Newton starts.  The first Jacobian callback
    // must therefore be a numeric update of existing slots, not structure discovery.
    PetscBool jacobianPreassembled = PETSC_FALSE;
    PetscCallAbort(
        mesh.communicator(),
        MatAssembled(jacobian, &jacobianPreassembled));
    if (!jacobianPreassembled)
    {
        throw std::runtime_error(
            "Natural Jacobian static topology was not materialized during matrix creation.");
    }

    SNES snes = nullptr;

    PetscCallAbort(
        mesh.communicator(),
        SNESCreate(
            mesh.communicator(),
            &snes));

    PetscCallAbort(
        mesh.communicator(),
        MPMC::installNaturalCallbacks(
            snes,
            runtime,
            residual,
            jacobian,
            jacobian));

    /*
     * 用户修改版 PETSc 的 updateState/updateSol 从 SNES application context
     * 取 Runtime，因此这里把这条上下文链也作为正式测试的一部分。
     */
    void *applicationContext = nullptr;
    PetscCallAbort(
        mesh.communicator(),
        SNESGetApplicationContext(
            snes,
            &applicationContext));

    if (applicationContext != &runtime)
    {
        throw std::runtime_error(
            "SNES application context is not the Natural runtime.");
    }

    const auto scatterBeforeFunction =
        runtime.vecScatterStatistics();

    // Normal residual assembly intentionally no longer performs VecSet(0): every
    // owned equation row is overwritten exactly once before well-control repair.
    // Seed NaN so any missed row is caught by the existing finite-vector check.
    PetscCallAbort(
        mesh.communicator(),
        VecSet(
            residual,
            static_cast<PetscScalar>(
                std::numeric_limits<PetscReal>::quiet_NaN())));

    PetscCallAbort(
        mesh.communicator(),
        SNESComputeFunction(
            snes,
            solution,
            residual));

    const auto scatterAfterFunction =
        runtime.vecScatterStatistics();
    if (scatterAfterFunction.currentSolution !=
            scatterBeforeFunction.currentSolution + 1 ||
        scatterAfterFunction.currentPhaseState !=
            scatterBeforeFunction.currentPhaseState + 1 ||
        scatterAfterFunction.maximumGasSaturation !=
            scatterBeforeFunction.maximumGasSaturation ||
        scatterAfterFunction.total() != scatterBeforeFunction.total() + 2)
    {
        throw std::runtime_error(
            "Natural CellCache miss did not use exactly solution+phase ghost scatters.");
    }

    PetscCallAbort(
        mesh.communicator(),
        SNESComputeJacobian(
            snes,
            solution,
            jacobian,
            jacobian));

    const auto scatterAfterJacobian =
        runtime.vecScatterStatistics();
    if (scatterAfterJacobian.total() != scatterAfterFunction.total())
    {
        throw std::runtime_error(
            "Natural Jacobian did not reuse the residual CellCache snapshot.");
    }

    // This integration model has Land trapping disabled.  A matching accepted
    // CellCache therefore lets commit reuse state/history with zero extra G2L.
    runtime.commitTimeStep(solution);
    const auto scatterAfterCommit =
        runtime.vecScatterStatistics();
    if (scatterAfterCommit.maximumGasSaturation !=
            scatterAfterJacobian.maximumGasSaturation ||
        scatterAfterCommit.previousSolution !=
            scatterAfterJacobian.previousSolution ||
        scatterAfterCommit.previousPhaseState !=
            scatterAfterJacobian.previousPhaseState ||
        scatterAfterCommit.total() != scatterAfterJacobian.total())
    {
        throw std::runtime_error(
            "Land-off accepted Natural commit performed an unnecessary ghost scatter.");
    }

    // 数值：trial-state EOS 异常必须成为可恢复的 SNES domain failure，
    // 不能以 PETSc error code 直接绕过自适应时间步重试。
    SNES domainSnes = nullptr;
    Vec domainResidual = nullptr;
    DomainFailureRuntime domainRuntime;
    PetscCallAbort(mesh.communicator(), SNESCreate(mesh.communicator(), &domainSnes));
    PetscCallAbort(mesh.communicator(), VecDuplicate(residual, &domainResidual));
    PetscCallAbort(
        mesh.communicator(),
        MPMC::naturalFormFunction<DomainFailureRuntime>(
            domainSnes, solution, domainResidual, &domainRuntime));
    PetscBool domainError = PETSC_FALSE;
    PetscCallAbort(
        mesh.communicator(),
        SNESGetFunctionDomainError(domainSnes, &domainError));
    if (domainError != PETSC_TRUE)
        throw std::runtime_error("Natural residual exception was not reported as a SNES domain error.");
    PetscReal domainResidualNorm = 0.0;
    PetscCallAbort(mesh.communicator(), VecNorm(domainResidual, NORM_INFINITY, &domainResidualNorm));
    if (domainResidualNorm != 0.0)
        throw std::runtime_error("Recoverable domain failure did not clear the residual vector.");
    PetscCallAbort(mesh.communicator(), VecDestroy(&domainResidual));
    PetscCallAbort(mesh.communicator(), SNESDestroy(&domainSnes));

    requireFiniteVector(
        residual,
        "Residual");

    MatInfo localInfo{};

    PetscCallAbort(
        mesh.communicator(),
        MatGetInfo(
            jacobian,
            MAT_LOCAL,
            &localInfo));

    PetscLogDouble globalMallocs = 0.0;
    const PetscLogDouble localMallocs =
        localInfo.mallocs;

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Allreduce(
            &localMallocs,
            &globalMallocs,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            mesh.communicator()));

    if (globalMallocs != 0.0)
    {
        throw std::runtime_error(
            "Natural Jacobian performed unexpected dynamic matrix allocation.");
    }

    PetscInt localCrossRankFaces = 0;

    for (const auto &cell :
         grid.localCells())
    {
        for (const auto &face :
             cell.faces())
        {
            const auto *neighbor =
                face.neighborCell();

            if (neighbor != nullptr &&
                neighbor->processorId() !=
                    cell.processorId())
            {
                ++localCrossRankFaces;
            }
        }
    }

    PetscInt globalCrossRankFaces = 0;

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Allreduce(
            &localCrossRankFaces,
            &globalCrossRankFaces,
            1,
            MPIU_INT,
            MPI_SUM,
            mesh.communicator()));

    if (mesh.processCount() > 1 &&
        globalCrossRankFaces == 0)
    {
        throw std::runtime_error(
            "Multi-rank test did not contain any cross-rank face.");
    }

    if (mesh.rank() == 0 && mesh.processCount() > 1)
    {
        const auto &representative = mesh.cellByCurrentId(0);
        const auto &remotePerforation = mesh.cellByCurrentId(
            static_cast<PetscInt>(mesh.cellCount() - 1));
        if (representative.processorId() == remotePerforation.processorId())
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "[WARN] first/last well perforations are on the same rank; "
                "topological cross-rank faces are still being tested.\n");
        }
    }

    const auto diagnostics = runtime.evaluateGlobalDiagnostics(solution);
    if (diagnostics.cellCount != static_cast<long long>(mesh.cellCount()) ||
        !std::isfinite(diagnostics.pressureAverage) ||
        !std::isfinite(diagnostics.porosityAverage))
    {
        throw std::runtime_error(
            "Natural/CpGrid global reservoir diagnostics are inconsistent.");
    }
    const double averageSaturationSum =
        diagnostics.saturationAverage[Indices::Phase::liquid] +
        diagnostics.saturationAverage[Indices::Phase::vapor] +
        diagnostics.saturationAverage[Indices::Phase::water];
    if (std::abs(averageSaturationSum - 1.0) > 1.0e-10)
        throw std::runtime_error("Natural/CpGrid average phase saturations do not sum to one.");

    const double jacobianRelativeError =
        checkPressureJacobianColumn(
            snes,
            runtime,
            solution,
            jacobian);

    /*
     * PR/LBC + upwind 模型的中心差分受浮点尺度影响，1e-4 对完整残差 Jacobian
     * 已足够严格，同时不会把非光滑切换点误判成 AD 错误。
     */
    if (!(jacobianRelativeError < 1.0e-4) ||
        !std::isfinite(jacobianRelativeError))
    {
        throw std::runtime_error(
            "AD/PETSc Jacobian finite-difference check failed: relative error = " +
            std::to_string(jacobianRelativeError));
    }

    /*
     * 单独实际调用 updateSol adapter。零 step 不改变解，但完整覆盖用户修改版
     * PETSc 所要求的符号和 Vec array 生命周期。
     */
    Vec newtonStep = nullptr;

    PetscCallAbort(
        mesh.communicator(),
        VecDuplicate(
            solution,
            &newtonStep));
    PetscCallAbort(
        mesh.communicator(),
        VecSet(
            newtonStep,
            0.0));

    runtime.limitNewtonStep(
        solution,
        newtonStep);

    requireFiniteVector(
        newtonStep,
        "Limited Newton step");

    PetscPrintf(
        mesh.communicator(),
        "Natural/CpGrid PETSc adapter validation\n"
        "  ranks                 = %d\n"
        "  cells                 = %zu\n"
        "  primary vars/cell     = %d\n"
        "  cross-rank face sides = %" PetscInt_FMT "\n"
        "  matrix mallocs        = %.0f\n"
        "  FD/AD Jacobian error  = %.3e\n"
        "  Vec G2L scatters      = %llu\n"
        "  result                = ALL PASS\n",
        mesh.processCount(),
        mesh.cellCount(),
        Indices::numPrimaryVariables,
        globalCrossRankFaces,
        static_cast<double>(
            globalMallocs),
        jacobianRelativeError,
        static_cast<unsigned long long>(
            runtime.vecScatterStatistics().total()));

    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(&newtonStep));
    PetscCallAbort(
        mesh.communicator(),
        SNESDestroy(&snes));
    PetscCallAbort(
        mesh.communicator(),
        MatDestroy(&jacobian));
    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(&residual));
    PetscCallAbort(
        mesh.communicator(),
        VecDestroy(&solution));
}

} // namespace AdapterTest

/**
 * @brief 用户修改版 PETSc (ls.c/virs.c/viss.c) 需要的 C++ linkage 外部符号。
 *
 * clean PETSc 不调用这两个符号；在用户修改版 PETSc 中，install helper
 * 通过 SNESSetApplicationContext() 显式把 Runtime 绑定到 snes->user 路径。
 */
void updateState(Vec solution, void *context)
{
    MPMC::naturalUpdateStateHook<AdapterTest::Runtime>(
        solution,
        context);
}

void updateSol(Vec solution, Vec step, void *context)
{
    MPMC::naturalUpdateSolHook<AdapterTest::Runtime>(
        solution,
        step,
        context);
}

int main(int argc, char **argv)
{
    PetscInitialize(
        &argc,
        &argv,
        nullptr,
        nullptr);

    int exitCode = 0;

    try
    {
        char meshDirectory[
            PETSC_MAX_PATH_LEN] = {};

        PetscBool found =
            PETSC_FALSE;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetString(
                nullptr,
                nullptr,
                "-mesh_dir",
                meshDirectory,
                sizeof(meshDirectory),
                &found));

        if (!found)
        {
            throw std::invalid_argument(
                "Please provide -mesh_dir <MRST-export-directory>.");
        }

        AdapterTest::runTest(
            std::string(
                meshDirectory));
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            PETSC_COMM_WORLD,
            "Natural/CpGrid PETSc adapter validation failed: %s\n",
            error.what());

        exitCode = 1;
    }

    PetscFinalize();
    return exitCode;
}
