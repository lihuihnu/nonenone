/**
 * @file natural_petsc_runtime.hpp
 * @brief 与网格后端解耦的 Natural/PETSc 非线性求解运行时。
 */
#pragma once

#include <common/math.hpp>
#include <natural/kernel/cell_kernel.hpp>
#include <natural/petsc/jacobian_factory.hpp>
#include <natural/petsc/grid_backend_common.hpp>
#include <natural/petsc/natural_runtime_options.hpp>
#include <natural/petsc/natural_runtime_diagnostics.hpp>
#include <natural/petsc/phase_state_codec.hpp>
#include <natural/petsc/property_conversion.hpp>
#include <natural/petsc/natural_layout_registration.hpp>
#include <natural/petsc/vector_access.hpp>
#include <natural/petsc/well_runtime.hpp>
#include <natural/physics/accumulation.hpp>
#include <natural/primary_variables.hpp>
#include <natural/state/newton_limiter.hpp>

#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace MPMC
{


/**
 * @brief 与具体网格无关的 Natural/PETSc 运行时。
 *
 * Backend 只封装 DOF/ghost/几何/岩石属性差异；EOS、flash、通量、蓄积、
 * 井源项、Newton limiter 和 Jacobian 装配流程只有这一份实现。
 */
template <
    class Indices,
    class Backend>
class NaturalPetscRuntime final
{
public:
    static_assert(
        Indices::numPrimaryVariables ==
            Indices::numEquations,
        "Natural PETSc runtime requires a square cell system.");

    using Scalar = typename Indices::ValueType;
    using Grid = typename Backend::GridType;
    using Cell = typename Backend::CellId;
    using DofMap = typename Backend::DofMap;
    using Kernel = NaturalCellKernel<Indices>;
    using State = CellState<Indices, Scalar>;
    using Properties = CellProperties<Indices, Scalar>;
    using ScalarProperties = CellProperties<Indices, double>;
    using PrimaryArray =
        std::array<double, Indices::numPrimaryVariables>;
    using Composition =
        std::array<Scalar, Indices::numComponents>;
    using Well = NaturalWell<Indices>;
    using WellStateType = WellState<Indices>;
    using Options = NaturalRuntimeOptions<Indices>;
    using ResidualFailureDiagnostic = NaturalMassBalanceFailureDiagnostic<Indices>;

    struct StateUpdateReport final
    {
        bool succeeded{true};
        std::size_t localUpdatedCellCount{0};
        std::size_t localStableReducedCellCount{0};
        std::size_t failingRankCount{0};
        int firstFailureRank{-1};
        Cell firstFailureCell{Cell(-1)};
        PhaseUpdateResult firstFailure{};
    };

    explicit NaturalPetscRuntime(
        Grid &grid,
        const FluidSystem<Indices> &fluid,
        Options options = {})
        : backend_(grid),
          fluid_(fluid),
          kernel_(fluid),
          options_(std::move(options))
    {
        validateOptions_();
        resolvedCO2ComponentCache_ = resolveCO2Component_();

        // Natural 显式拥有模型布局注册职责；Grid 核心只注册自身几何/岩石布局。
        registerNaturalGridLayouts<Indices>(backend_);

        phaseState_ =
            backend_.createGlobalVector(
                static_cast<PetscInt>(
                    Indices::numPhaseStateVariables));

        previousSolution_ =
            backend_.createGlobalVector(
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables));

        previousPhaseState_ =
            backend_.createGlobalVector(
                static_cast<PetscInt>(
                    Indices::numPhaseStateVariables));

        maximumGasSaturation_ =
            backend_.createGlobalVector(1);

        PetscCallAbort(
            backend_.communicator(),
            VecSet(phaseState_, 0.0));
        PetscCallAbort(
            backend_.communicator(),
            VecSet(previousSolution_, 0.0));
        PetscCallAbort(
            backend_.communicator(),
            VecSet(previousPhaseState_, 0.0));
        PetscCallAbort(
            backend_.communicator(),
            VecSet(maximumGasSaturation_, 0.0));

        initializeStaticAssemblyCache_();
    }

    NaturalPetscRuntime(
        const NaturalPetscRuntime &) = delete;
    NaturalPetscRuntime &operator=(
        const NaturalPetscRuntime &) = delete;
    NaturalPetscRuntime(
        NaturalPetscRuntime &&) = delete;
    NaturalPetscRuntime &operator=(
        NaturalPetscRuntime &&) = delete;

    ~NaturalPetscRuntime() noexcept
    {
        destroyVector_(phaseState_);
        destroyVector_(previousSolution_);
        destroyVector_(previousPhaseState_);
        destroyVector_(maximumGasSaturation_);
    }

    /** @brief 访问由调用方拥有的网格后端对象；Runtime 不接管其生命周期。 */
    [[nodiscard]] Grid &grid() noexcept
    {
        return backend_.grid();
    }

    /** @brief 只读访问由调用方拥有的网格后端对象。 */
    [[nodiscard]] const Grid &grid() const noexcept
    {
        return backend_.grid();
    }

    [[nodiscard]] PetscInt cellCount() const noexcept
    {
        return backend_.cellCount();
    }

    [[nodiscard]] const FluidSystem<Indices> &fluidSystem() const noexcept
    {
        return fluid_;
    }

    /** @brief 返回当前 rank 累计触发的 Natural global-to-local Vec 同步次数。 */
    [[nodiscard]] const NaturalVecScatterStatistics &vecScatterStatistics() const noexcept
    {
        return vecScatterStatistics_;
    }


    /**
     * @brief 创建按网格标准输入顺序排列的单元 Vec 副本。
     *
     * CpGrid solves with partition/current-id ordering internally, while external
     * files and well tables use the original grid-file row index.  StructuredGrid
     * has no extra permutation, so this operation is just a copy there.  The caller
     * owns the returned Vec and must destroy it.
     */
    [[nodiscard]] Vec createInputOrderedCopy(
        Vec currentOrderedVector,
        PetscInt dofPerCell) const
    {
        return backend_.createInputOrderedCopy(
            dofPerCell,
            currentOrderedVector);
    }

    /**
     * @brief 返回内部单元 id 对应的标准外部/输入单元序号。
     */
    [[nodiscard]] PetscInt inputCellIndex(Cell currentCellId) const
    {
        return backend_.inputCellIndex(currentCellId);
    }

#include <natural/petsc/detail/natural_petsc_runtime_failure_diagnostics.inc>

    /** @brief 保存 K、z、L、Z 等二级相态量的 PETSc 全局向量。 */
    [[nodiscard]] Vec phaseStateVector() const noexcept
    {
        return phaseState_;
    }

    /** @brief 当前时间步开始时已经接受的主变量解，用于时间离散和回滚。 */
    [[nodiscard]] Vec previousSolutionVector() const noexcept
    {
        return previousSolution_;
    }

    /** @brief 当前时间步开始时已经接受的相态，用于失败尝试回滚。 */
    [[nodiscard]] Vec previousPhaseStateVector() const noexcept
    {
        return previousPhaseState_;
    }

    /** @brief 设置蓄积项使用的全隐式时间步长。 */
    void setTimeStep(double value)
    {
        if (!(value > 0.0) ||
            !std::isfinite(value))
        {
            throw std::invalid_argument(
                "Natural time step must be finite and positive.");
        }

        options_.timeStep = value;
    }

    /** @brief 校验穿孔单元和 BHP 代表单元后，整体替换当前井集合。 */
    void setWells(std::vector<Well> wells)
    {
        if constexpr (!Indices::hasWellUnknown)
        {
            if (!wells.empty())
            {
                throw std::invalid_argument(
                    "This Indices configuration has no well unknown.");
            }
        }

        const PetscInt cellCount =
            static_cast<PetscInt>(
                backend_.cellCount());

        std::unordered_set<PetscInt>
            representativeCells;

        for (const auto &well : wells)
        {
            well.validate(cellCount);

            if (!representativeCells.insert(
                     well.bhpCellId)
                     .second)
            {
                throw std::invalid_argument(
                    "Two wells cannot share the same BHP representative cell.");
            }
        }

        cacheWellGlobalDofs_(wells);
        cacheWellAssemblyMetadata_(wells);
        wells_ = std::move(wells);
        jacobianPatternValid_ = false;
        cachedWellBhpCellCacheGeneration_ =
            std::numeric_limits<std::size_t>::max();
    }

    /** @brief 以稳定插入顺序只读访问井定义。 */
    [[nodiscard]] const std::vector<Well> &wells() const noexcept
    {
        return wells_;
    }

    /** @brief 可写访问井定义，仅供收敛后的井控制更新使用。 */
    [[nodiscard]] std::vector<Well> &mutableWells() noexcept
    {
        return wells_;
    }

    using WellScheduleUpdater =
        std::function<void(double, std::vector<Well> &)>;

    /**
     * @brief 安装可选的随时间更新井运行模式的回调。
     *
     * This hook is intentionally more general than WellSchedule's simple
     * open/close window.  It allows a physical well to change role/control
     * (for example producer -> idle -> injector) without duplicating the well
     * and therefore without allocating two BHP equations to the same cell.
     * The updater is evaluated whenever the runtime time changes, including
     * adaptive trial end times.
     */
    void setWellScheduleUpdater(WellScheduleUpdater updater)
    {
        wellScheduleUpdater_ = std::move(updater);
        if (wellScheduleUpdater_)
            wellScheduleUpdater_(currentTime_, wells_);
    }

    /** @brief 设置 schedule 与井启停判断使用的物理模拟时间。 */
    void setCurrentTime(double time)
    {
        if (!std::isfinite(time))
            throw std::invalid_argument(
                "Natural current time must be finite.");
        currentTime_ = time;
        if (wellScheduleUpdater_)
            wellScheduleUpdater_(currentTime_, wells_);
    }

    /** @brief 返回当前物理模拟时间。 */
    [[nodiscard]] double currentTime() const noexcept
    {
        return currentTime_;
    }

#include <natural/petsc/detail/natural_petsc_runtime_output_diagnostics.inc>

    [[nodiscard]] Vec createResidualVector() const
    {
        return backend_.createGlobalVector(
            static_cast<PetscInt>(
                Indices::numEquations));
    }

    /**
     * @brief 创建包含拓扑块和井非局部块的 Jacobian。
     */
    [[nodiscard]] Mat createJacobian() const
    {
        if (!jacobianPatternValid_)
        {
            jacobianPatternCache_ =
                buildNaturalJacobianPattern<Indices>(backend_, wells_);
            jacobianPatternValid_ = true;
        }

        return createNaturalJacobianFromPattern(
            backend_,
            jacobianPatternCache_);
    }

#include <natural/petsc/detail/natural_petsc_runtime_state.inc>

#include <natural/petsc/detail/natural_petsc_runtime_residual.inc>

private:
    struct CellCacheEntry final
    {
        State state{};
        Composition overallComposition{};
        Properties properties{};
    };

    /**
     * @brief residual/Jacobian 热循环直接消费的静态连接。
     *
     * neighborLocalBlock 在 runtime 构造时解析一次，避免每次面通量都把
     * neighbor current-id 再映射到 local Vec block。
     */
    struct StaticAssemblyConnection final
    {
        Cell neighborCellId{Cell(-1)};
        PetscInt neighborLocalBlock{-1};
        double transmissibility{0.0};
        double gravityTerm{0.0};
    };

    /** @brief 当前 rank 真正拥有的一条井穿孔的静态装配索引。 */
    struct LocalWellPerforationAssembly final
    {
        std::size_t perforationIndex{0};
        PetscInt localBlock{-1};
        PetscInt globalBlock{-1};
    };

    /** @brief 每口井在当前 rank 上不随 Newton 改变的装配元数据。 */
    struct WellStaticAssemblyCache final
    {
        PetscInt representativeLocalBlock{-1};
        std::vector<LocalWellPerforationAssembly> localPerforations;
    };

    struct CellWellSource final
    {
        std::array<Scalar, Indices::numComponents>
            component{};
        Scalar water{0.0};
    };

    struct WellFunctionContributions final
    {
        std::vector<CellWellSource> cellSource;
        std::vector<double> controlledRateGlobal;
    };

    struct PendingStateUpdate final
    {
        Cell cell{Cell(-1)};
        PrimaryArray primary{};
        PhaseStateData<Indices> phaseState{};
    };

#include <natural/petsc/detail/natural_petsc_runtime_common.inc>

#include <natural/petsc/detail/natural_petsc_runtime_cell_cache.inc>

#include <natural/petsc/detail/natural_petsc_runtime_wells.inc>

#include <natural/petsc/detail/natural_petsc_runtime_jacobian.inc>

#include <natural/petsc/detail/natural_petsc_runtime_scaling.inc>

#include <natural/petsc/detail/natural_petsc_runtime_cache_history.inc>

    Backend backend_;
    const FluidSystem<Indices> &fluid_;
    Kernel kernel_;
    Options options_;
    double currentTime_{0.0};
    WellScheduleUpdater wellScheduleUpdater_{};

    std::vector<Well> wells_;
    std::vector<PetscInt> wellBhpGlobalDof_;
    std::vector<PetscInt> wellControlGlobalDof_;
    std::vector<WellStaticAssemblyCache> wellAssemblyCache_;
    std::vector<PetscInt> localWellSourceBlocks_;
    std::vector<std::size_t> representativeWellByLocalBlock_;

    std::vector<Cell> ownedCells_;
    std::vector<PetscInt> ownedLocalBlockCache_;
    std::vector<Cell> localSnapshotCells_;
    std::vector<PetscInt> localSnapshotBlockCache_;
    std::size_t localCacheBlockCount_{0};
    std::vector<double> porosityCache_;
    std::vector<double> cellVolumeCache_;
    std::vector<double> maximumGasLocalCache_;
    std::vector<std::vector<StaticAssemblyConnection>> connectionCache_;
    // G8K: residual 只保存每条有向连接的各相 upwind 决策，不缓存逐面 AD。
    // generation 按 owned local block 记录；Jacobian 只在状态完全相同时复用。
    mutable std::vector<std::vector<FaceFluxLinearizationDecision<Indices>>>
        faceFluxDecisionCache_;
    mutable std::vector<std::size_t> faceFluxDecisionGeneration_;
    std::vector<PetscInt> primaryGlobalBlockCache_;
    std::vector<ScalarProperties> previousPropertiesCache_;
    std::vector<AccumulationResult<Indices, double>> previousFluidAccumulationCache_;
    std::vector<std::array<double, Indices::numComponents>> previousAdsorbedAccumulationCache_;
    mutable std::vector<CellCacheEntry> cellCacheScratch_;
    mutable std::vector<double> wellBhpLocalScratch_;
    mutable std::vector<double> wellBhpGlobalScratch_;
    mutable std::vector<Scalar> wellBhpScalarScratch_;
    mutable WellFunctionContributions wellFunctionScratch_;
    mutable std::vector<double> wellRateLocalScratch_;
    mutable std::array<Scalar, Indices::numPhases> wellPhasePressureScratch_{};
    mutable PerforationWellResult<Indices, Scalar> wellPerforationResultScratch_{};
    mutable PerforationWellWorkspace<Indices, Scalar> wellPerforationWorkspaceScratch_{};
    // 一个 owned cell 的 Jacobian 对同一 interior column block 写入：
    // [self row, reciprocal neighbor rows...]。scratch 在 setup 时按最大连接数
    // 一次性定长，Newton 热循环只覆盖有效前缀，不再 clear/resize。
    mutable std::vector<PetscInt> cellJacobianRowBlockScratch_;
    mutable std::vector<PetscScalar> cellJacobianBlockValueScratch_;
    mutable NaturalJacobianPattern jacobianPatternCache_{};
    mutable bool jacobianPatternValid_{false};
    mutable PetscInt preparedCellJacobianRowCount_{0};
    mutable std::size_t cellCacheGeneration_{0};
    mutable std::size_t cachedWellBhpCellCacheGeneration_{
        std::numeric_limits<std::size_t>::max()};
    mutable Vec cachedSolution_{nullptr};
    mutable PetscObjectState cachedSolutionState_{};
    mutable PetscObjectState cachedPhaseStateState_{};
    mutable PetscObjectState cachedMaximumGasState_{};
    mutable bool currentCacheValid_{false};
    mutable NaturalVecScatterStatistics vecScatterStatistics_{};
    // PETSc < 3.24 compatibility: synthetic ever-changing Vec state disables
    // cache reuse when no public Vec state-query API is available.
    mutable PetscObjectState compatibilityVectorStateCounter_{};

    int resolvedCO2ComponentCache_{-1};

    Vec phaseEquilibriumSolution_{nullptr};
    PetscObjectState phaseEquilibriumSolutionState_{};
    PetscObjectState phaseEquilibriumPhaseStateState_{};
    bool phaseEquilibriumCacheValid_{false};

    std::vector<PendingStateUpdate> stateUpdateScratch_;
    StateUpdateReport lastStateUpdateReport_{};
    bool pendingThermoStateFailure_{false};

    Vec phaseState_{nullptr};
    Vec previousSolution_{nullptr};
    Vec previousPhaseState_{nullptr};
    Vec maximumGasSaturation_{nullptr};

    bool phaseStateInitialized_{false};
    bool historyInitialized_{false};
};

} // namespace MPMC
