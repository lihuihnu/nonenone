/**
 * @file cell_kernel.hpp
 * @brief 单元状态更新、物性评价和局部残差计算的统一计算内核。
 */
#pragma once

#include <natural/assembly/cell_residual.hpp>
#include <natural/fluid_system.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/phase_equilibrium.hpp>
#include <natural/state/three_phase_equilibrium.hpp>
#include <natural/state/state_codec.hpp>

#include <array>
#include <type_traits>

namespace MPMC
{

/**
 * @brief natural formulation 的单控制体纯物理入口。
 *
 * 这是上层 Simulator/PETSc 适配器推荐依赖的对象。它把原来散落在
 * IntensiveQuantities、FormFunction、updateState 中的纯局部物理统一到一个
 * 无网格/无 PETSc 的 API：
 *
 * `primary -> CellState -> CellProperties -> CellResidualBlock`。
 *
 * 分布式存储、ghost 更新、Mat/Vec 装配仍由外层 adapter 负责。
 */
template <class Indices>
class NaturalCellKernel final
{
public:
    using Scalar = typename Indices::ValueType;
    using PrimaryArray =
        std::array<double, Indices::numPrimaryVariables>;
    using State =
        CellState<Indices, Scalar>;
    using Properties =
        CellProperties<Indices, Scalar>;
    using ScalarProperties =
        CellProperties<Indices, double>;
    using Residual =
        CellResidualBlock<Indices, Scalar>;
    using Composition =
        std::array<Scalar, Indices::numComponents>;
    using PhaseEquilibrium = std::conditional_t<
        Indices::fullyCompositionalThreePhase,
        FullyCompositionalThreePhaseEquilibrium<Indices>,
        PhaseEquilibriumManager<Indices>>;

    explicit NaturalCellKernel(
        const FluidSystem<Indices> &fluid)
        : fluid_(fluid),
          propertyEvaluator_(fluid),
          phaseEquilibrium_(fluid)
    {
    }

    /** @brief 将 PETSc/算例侧主变量和二级相态解码为强类型单元状态。 */
    [[nodiscard]] State decodeState(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState) const
    {
        return CellStateCodec<Indices>::decode(
            primary,
            phaseState);
    }

    /** @brief 由单元状态统一评价 EOS、相物性、孔隙度修正和可选滞留/吸附量。 */
    [[nodiscard]] Properties evaluateProperties(
        const State &state,
        const Scalar &basePorosity,
        double maximumGasSaturation = 0.0,
        double landConstant = 2.0,
        const Composition *overallComposition = nullptr) const
    {
        return propertyEvaluator_.evaluate(
            state,
            basePorosity,
            maximumGasSaturation,
            landConstant,
            overallComposition);
    }

    /** @brief 组装仅依赖当前单元与上一时间层的 local/base residual。 */
    [[nodiscard]] Residual assembleLocalResidual(
        const State &state,
        const Properties &properties,
        const ScalarProperties &previousProperties,
        double cellVolume,
        double timeStep,
        double fugacityScalingFactor,
        int dissolvedCO2Component = -1,
        double rockDensity = 2650.0,
        const std::array<double, Indices::numComponents> *standardGasDensity = nullptr,
        const AccumulationResult<Indices, double> *previousFluidAccumulation = nullptr,
        const std::array<double, Indices::numComponents> *previousAdsorbedAccumulation = nullptr) const
    {
        return assembleCellLocalResidual<Indices, Scalar>(
            state,
            properties,
            previousProperties,
            cellVolume,
            timeStep,
            fugacityScalingFactor,
            dissolvedCO2Component,
            rockDensity,
            standardGasDensity,
            previousFluidAccumulation,
            previousAdsorbedAccumulation);
    }

    /** @brief 按统一守恒/相平衡离散组装单元完整残差；不负责网格或 PETSc 写入。 */
    [[nodiscard]] Residual assembleResidual(
        const State &state,
        const Properties &properties,
        const ScalarProperties &previousProperties,
        const FaceMassFlux<Indices, Scalar> &outwardFaceFlux,
        const std::array<Scalar, Indices::numComponents> &wellComponentMassSource,
        const Scalar &wellWaterMassSource,
        double cellVolume,
        double timeStep,
        double fugacityScalingFactor,
        int dissolvedCO2Component = -1,
        double rockDensity = 2650.0,
        const std::array<double, Indices::numComponents> *standardGasDensity = nullptr) const
    {
        return assembleCellResidual<Indices, Scalar>(
            state,
            properties,
            previousProperties,
            outwardFaceFlux,
            wellComponentMassSource,
            wellWaterMassSource,
            cellVolume,
            timeStep,
            fugacityScalingFactor,
            dissolvedCO2Component,
            rockDensity,
            standardGasDensity);
    }

    /** @brief 访问相态管理器；相消失/再出现的状态转换由该对象统一处理。 */
    [[nodiscard]] const PhaseEquilibrium &
    phaseEquilibrium() const noexcept
    {
        return phaseEquilibrium_;
    }

private:
    const FluidSystem<Indices> &fluid_;
    CellPropertyEvaluator<Indices> propertyEvaluator_;
    PhaseEquilibrium phaseEquilibrium_;
};

} // namespace MPMC
