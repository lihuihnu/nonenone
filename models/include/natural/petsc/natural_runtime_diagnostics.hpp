/**
 * @file natural_runtime_diagnostics.hpp
 * @brief Natural/PETSc runtime 的诊断与全局统计数据结构。
 */
#pragma once

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace MPMC
{


/**
 * @brief Natural runtime 实际触发的 global-to-local Vec 同步计数。
 *
 * 这些计数只在 runtime 明确创建 local read snapshot 时递增，因此可用于
 * 区分“CellCache 命中”与真正发生的 PETSc ghost scatter。计数是每个 MPI rank
 * 的本地累计值，不做隐式 MPI 归约。
 */
struct NaturalVecScatterStatistics final
{
    std::uint64_t currentSolution{0};
    std::uint64_t currentPhaseState{0};
    std::uint64_t previousSolution{0};
    std::uint64_t previousPhaseState{0};
    std::uint64_t maximumGasSaturation{0};

    [[nodiscard]] std::uint64_t total() const noexcept
    {
        return currentSolution +
               currentPhaseState +
               previousSolution +
               previousPhaseState +
               maximumGasSaturation;
    }
};

/**
 * @brief MPI 全局井组分源汇率。
 *
 * 注入量与产出量均保存为非负工程量（kg/s），netSource 则遵循储层
 * 残差中的符号约定：注入为正、产出为负。
 */
template <class Indices>
struct NaturalGlobalWellComponentRates final
{
    // 单位均为 kg/s。injected/produced 保存非负工程量；netSource 使用
    // 储层方程符号约定（注入为正、产出为负）。
    std::array<double, Indices::numComponents> netSource{};
    std::array<double, Indices::numComponents> injected{};
    std::array<double, Indices::numComponents> produced{};

    // 旧模型把 H2O 作为独立守恒方程；全组分 O/G/W 模型已把 H2O 放入
    // numComponents，因此以下独立水统计保持为零。
    double waterNetSource{0.0};
    double waterInjected{0.0};
    double waterProduced{0.0};
};

/**
 * @brief MPI 全局库存统计。
 *
 * 所有数值单位均为 kg。trappedGasMass 是气相总质量中被 Land 模型
 * 判定为残余滞留的那一部分，不应再额外加到 fluidComponentMass 中。
 */
template <class Indices>
struct NaturalGlobalInventory final
{
    std::array<double, Indices::numComponents> fluidComponentMass{};
    std::array<double, Indices::numComponents> oilPhaseComponentMass{};
    double waterMass{0.0};

    double dissolvedCO2Mass{0.0};

    double trappedGasMass{0.0};
    std::array<double, Indices::numComponents> trappedComponentMass{};

    std::array<double, Indices::numComponents> adsorbedComponentMass{};
};

/**
 * @brief 收敛状态的 MPI 全局储层诊断量。
 *
 * 平均值使用与物理量一致的权重：压力/孔隙度按 bulk volume，饱和度按
 * pore volume，相密度/黏度按对应相 pore volume。密度和黏度的极值统计
 * 忽略单元内饱和度接近零的相，避免消失相物性污染统计。
 */
template <class Indices>
struct NaturalGlobalDiagnostics final
{
    long long cellCount{0};
    double bulkVolume{0.0};
    double poreVolume{0.0};

    double pressureMinimum{0.0};
    double pressureAverage{0.0};
    double pressureMaximum{0.0};
    double porosityMinimum{0.0};
    double porosityAverage{0.0};
    double porosityMaximum{0.0};

    std::array<double, Indices::numPhases> phasePoreVolume{};
    std::array<double, Indices::numPhases> saturationMinimum{};
    std::array<double, Indices::numPhases> saturationAverage{};
    std::array<double, Indices::numPhases> saturationMaximum{};
    std::array<double, Indices::numPhases> densityMinimum{};
    std::array<double, Indices::numPhases> densityAverage{};
    std::array<double, Indices::numPhases> densityMaximum{};
    std::array<double, Indices::numPhases> viscosityMinimum{};
    std::array<double, Indices::numPhases> viscosityAverage{};
    std::array<double, Indices::numPhases> viscosityMaximum{};

    // 旧油/气相态计数，保持字段不变以兼容既有输出。
    long long liquidOnlyCells{0};
    long long vaporOnlyCells{0};
    long long twoPhaseCells{0};

    // 全组分相态按 PhasePresence::bits() 的 1..7 编码计数；索引 0 保留不用。
    std::array<long long, 8> phasePresenceCells{};

    // 由相边界滞回暂时抑制的缺失相单元数。它只反映 active-set 历史，
    // 不改变公开 phasePresence 编码；用于确认相消失后没有在边界立即 chatter。
    std::array<long long, 3> phaseSuppressionCells{};

    double aqueousCO2MassFractionAverage{0.0};
    double aqueousCO2MassFractionMaximum{0.0};
    double trappedGasSaturationAverage{0.0};
    double trappedGasSaturationMaximum{0.0};
};


/**
 * @brief 失败 Newton 的完整 scaled residual 最大行诊断。
 *
 * residualNorm2 / residualNormInfinity / maximumAbsoluteScaledResidual 都是
 * PETSc 实际收敛判据看到的“缩放后”残差。residualRms 用全局方程数归一化，
 * 用于判断固定 L2 absolute tolerance 是否存在网格尺寸效应。
 *
 * winning row 同时保存 equationScale 和反缩放后的 residual，便于区分
 * mass/fugacity/closure/well-control 各方程自身量纲。该结构只做诊断。
 */
template <class Indices>
struct NaturalFullResidualFailureDiagnostic final
{
    bool valid{false};
    long long globalEquationCount{0};
    double residualNorm2{0.0};
    double residualRms{0.0};
    double residualNormInfinity{0.0};

    double maximumAbsoluteScaledResidual{0.0};
    double signedScaledResidual{0.0};
    double equationScale{1.0};
    double signedUnscaledResidual{0.0};

    std::array<double, Indices::numComponents>
        globalSignedComponentMassResidual{};
    double globalSignedIndependentWaterResidual{0.0};

    PetscInt currentCellId{-1};
    PetscInt inputCellId{-1};
    int equationIndex{-1};
    double pressure{0.0};
    std::array<double, Indices::numPhases> saturation{};
    std::array<std::array<double, Indices::numComponents>, Indices::numPhases>
        moleFraction{};
    std::uint8_t phasePresenceBits{0};
    std::uint8_t phaseSuppressionBits{0};
};


/**
 * @brief MPI 全局有符号守恒残差 [kg/s]。
 *
 * component[] 对应 EOS/全组分质量方程。legacy independent-water 模型额外
 * 使用 independentWater；全组分 O/G/W 中 H2O 已经包含在 component[]。
 */
template <class Indices>
struct NaturalGlobalSignedMassResidual final
{
    std::array<double, Indices::numComponents> component{};
    double independentWater{0.0};

    [[nodiscard]] double maximumAbsolute() const noexcept
    {
        double maximum = 0.0;
        for (double value : component)
            maximum = std::max(maximum, std::abs(value));
        if constexpr (Indices::hasIndependentWaterConservation)
            maximum = std::max(maximum, std::abs(independentWater));
        return maximum;
    }
};


/**
 * @brief 失败 Newton 中一个面的守恒/上游诊断。
 *
 * 所有通量均采用“从诊断单元流出为正”的残差符号约定。
 */
template <class Indices>
struct NaturalMassBalanceFaceDiagnostic final
{
    PetscInt neighborCellId{-1};
    PetscInt neighborInputCellId{-1};
    double transmissibility{0.0};
    double gravityTerm{0.0};
    double neighborPressure{0.0};
    std::array<double, Indices::numPhases> neighborSaturation{};

    std::array<double, Indices::numPhases> potentialDifference{};
    // 0 = interior/diagnostic cell, 1 = exterior/neighbor cell.
    std::array<int, Indices::numPhases> upwindSide{};
    std::array<double, Indices::numPhases> upwindMobility{};
    std::array<double, Indices::numPhases> upwindDensity{};
    std::array<double, Indices::numPhases> darcyVolumeRate{};
    std::array<double, Indices::numPhases> phaseMassRate{};
    std::array<std::array<double, Indices::numComponents>, Indices::numPhases>
        phaseComponentMassFlux{};
    std::array<double, Indices::numComponents> totalComponentMassFlux{};
};

/**
 * @brief 失败 Newton 中诊断单元上的单个井穿孔贡献。
 */
template <class Indices>
struct NaturalMassBalanceWellDiagnostic final
{
    int wellVectorIndex{-1};
    int wellId{-1};
    int wellType{0};
    int wellControl{0};
    int active{0};
    double target{0.0};
    double signedTarget{0.0};
    double controlledRate{0.0};
    double controlResidual{0.0};
    double bhp{0.0};
    double wellIndex{0.0};
    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    std::array<double, Indices::numComponents> injectionComponentMassFraction{};
    std::array<double, Indices::numPhases> pressureDrop{};
    std::array<double, Indices::numPhases> mobility{};
    std::array<double, Indices::numPhases> density{};
    std::array<double, Indices::numPhases> surfacePhaseRate{};
    std::array<double, Indices::numPhases> reservoirPhaseRate{};
    std::array<double, Indices::numPhases> phaseMassRate{};
    std::array<double, Indices::numComponents> componentMassSource{};
};

/**
 * @brief 失败 Newton 的组分质量守恒分解诊断。
 *
 * 诊断对象是 MPI 全局绝对值最大的“组分质量守恒”残差，而不是相平衡
 * 方程。对该单元重新按生产 residual 路径分解
 * `R_i = accumulation_i + sum(face_flux_i) - well_source_i`，并保留每个
 * 邻接面每一相的 potential/upwind/mobility/flux 以及每个井穿孔的源项。
 * 该结构只用于失败诊断，不参与任何求解路径。
 */
template <class Indices>
struct NaturalMassBalanceFailureDiagnostic final
{
    bool valid{false};
    double maximumAbsoluteMassResidual{0.0};
    PetscInt currentCellId{-1};
    PetscInt inputCellId{-1};
    int componentIndex{-1};
    double timeStep{0.0};
    double cellVolume{0.0};
    double pressure{0.0};
    std::array<double, Indices::numPhases> saturation{};
    std::array<double, Indices::numPhases> density{};
    std::array<double, Indices::numPhases> mobility{};
    std::array<std::array<double, Indices::numComponents>, Indices::numPhases>
        massFraction{};

    // kg/m^3 bulk volume before multiplication by V/dt.
    std::array<double, Indices::numComponents> currentAccumulationDensity{};
    std::array<double, Indices::numComponents> previousAccumulationDensity{};

    // kg/s in the residual sign convention.
    std::array<double, Indices::numComponents> accumulationTerm{};
    std::array<double, Indices::numComponents> faceFluxTerm{};
    std::array<double, Indices::numComponents> wellSourceTerm{};
    std::array<double, Indices::numComponents> reconstructedResidual{};
    std::array<double, Indices::numComponents> petscResidual{};
    std::array<double, Indices::numComponents> reconstructionError{};

    std::vector<NaturalMassBalanceFaceDiagnostic<Indices>> faces;
    std::vector<NaturalMassBalanceWellDiagnostic<Indices>> perforations;
};


} // namespace MPMC
