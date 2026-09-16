/**
 * @file cell_residual.hpp
 * @brief Natural 单元守恒方程残差的局部装配。
 */
#pragma once

#include <common/math.hpp>
#include <common/units.hpp>
#include <natural/assembly/local_equations.hpp>
#include <natural/numerics.hpp>
#include <natural/physics/accumulation.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/state/cell_state.hpp>

#include <array>
#include <cstddef>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 单控制体完整 residual block。
 *
 * `value[equation]` 的排列严格由 Indices::Equation 决定。
 */
template <class Indices, class Scalar>
struct CellResidualBlock
{
    std::array<Scalar, Indices::numEquations> value{};
};

/**
 * @brief 清除 RATE 井代表单元的默认 well-control 占位行。
 *
 * Natural 在每个单元都保留一个 wellPressure 主变量。普通单元需要
 * `wellPressure - value(wellPressure)=0` 的单位 Jacobian 来闭合这个 dummy
 * 未知量；RATE 井的 BHP 代表单元则会把该 residual 行替换为 `q-target`，
 * 因而 Jacobian 也必须同步替换，不能保留默认的 +1 对角。
 */
template <class Indices, class Scalar>
void clearRateWellControlPlaceholder(CellResidualBlock<Indices, Scalar> &block)
{
    static_assert(Indices::hasWellUnknown,
                  "RATE well-control placeholder requires a well-pressure unknown.");
    block.value[static_cast<std::size_t>(Indices::Equation::wellControl)] = Scalar(0.0);
}

/**
 * @brief 组装单个控制体的 local/base 方程，不包含 PETSc、面通量或井源。
 *
 * 守恒行在这里仅包含 `(M^{n+1}-M^n)V_b/dt` 时间项；跨单元面通量和
 * 井源由 `addCellConservationTransportAndWell()` 在外层追加。吸附质量与
 * 流体质量使用同一个时间离散。
 *
 * 方程包括：
 * 1. 组分质量守恒；
 * 2. 水质量守恒；
 * 3. 井变量占位方程（零残差 + BHP identity Jacobian）；
 * 4. 两相逸度平衡或单相退化方程；
 * 5. 饱和度闭合；
 * 6. 可选水相 CO2 逸度平衡。
 *
 * @param currentState 当前 Newton 状态（可为 AD）。
 * @param currentProperties 当前局部物性（可为 AD）。
 * @param previousProperties 上一时间层物性，始终视为常量。
 * @param cellVolume 控制体体积。
 * @param timeStep 时间步长。
 * @param fugacityScalingFactor 逸度方程的无量纲缩放因子。
 */
template <class Indices,
          class Scalar>
[[nodiscard]] CellResidualBlock<Indices, Scalar>
assembleCellLocalResidual(
    const CellState<Indices, Scalar> &currentState,
    const CellProperties<Indices, Scalar> &currentProperties,
    const CellProperties<Indices, double> &previousProperties,
    double cellVolume,
    double timeStep,
    double fugacityScalingFactor,
    int dissolvedCO2Component = -1,
    [[maybe_unused]] double rockDensity = 2650.0,
    [[maybe_unused]] const std::array<double, Indices::numComponents> *standardGasDensity = nullptr,
    const AccumulationResult<Indices, double> *previousFluidAccumulation = nullptr,
    const std::array<double, Indices::numComponents> *previousAdsorbedAccumulation = nullptr)
{
    if (!(cellVolume > 0.0))
        throw std::invalid_argument(
            "Cell volume must be positive.");

    if (!(timeStep > 0.0))
        throw std::invalid_argument(
            "Time step must be positive.");

    if (!(fugacityScalingFactor > 0.0))
        throw std::invalid_argument(
            "Fugacity scaling factor must be positive.");

    CellResidualBlock<Indices, Scalar> residual;

    const double storageFactor =
        cellVolume / timeStep;

    AccumulationResult<Indices, double> previousAccumulationStorage{};
    if (previousFluidAccumulation == nullptr)
    {
        previousAccumulationStorage =
            computeFluidAccumulation<Indices>(
                previousProperties,
                dissolvedCO2Component);
        previousFluidAccumulation = &previousAccumulationStorage;
    }

    for (int component = 0;
         component < Indices::numComponents;
         ++component)
    {
        const std::size_t c =
            static_cast<std::size_t>(component);

        const int equation =
            Indices::Equation::massConservation[c];

        residual.value[
            static_cast<std::size_t>(equation)] =
            (computeFluidComponentAccumulation<Indices>(
                 currentProperties, component, dissolvedCO2Component) -
             previousFluidAccumulation->componentMass[c]) *
                storageFactor;
    }

    if constexpr (Indices::hasAdsorption)
    {
        if (standardGasDensity == nullptr)
            throw std::invalid_argument(
                "Adsorption requires standard gas densities.");

        std::array<double, Indices::numComponents> previousAdsorptionStorage{};
        if (previousAdsorbedAccumulation == nullptr)
        {
            previousAdsorptionStorage =
                computeAdsorbedAccumulation<Indices>(
                    previousProperties.porosity,
                    rockDensity,
                    *standardGasDensity,
                    previousProperties.adsorbedVolume);
            previousAdsorbedAccumulation = &previousAdsorptionStorage;
        }

        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            const std::size_t c =
                static_cast<std::size_t>(component);

            const int equation =
                Indices::Equation::massConservation[c];

            residual.value[
                static_cast<std::size_t>(equation)] +=
                (computeAdsorbedComponentAccumulation<Indices>(
                     currentProperties.porosity,
                     rockDensity,
                     *standardGasDensity,
                     currentProperties.adsorbedVolume,
                     component) -
                 (*previousAdsorbedAccumulation)[c]) *
                storageFactor;
        }
    }
    if constexpr (Indices::hasIndependentWaterConservation)
    {
        residual.value[
            static_cast<std::size_t>(
                Indices::Equation::waterConservation)] =
            (computeFluidWaterAccumulation<Indices>(currentProperties) -
             previousFluidAccumulation->waterMass) *
                storageFactor;
    }

    // 相平衡/相消失方程块。
    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        const bool oil = currentState.phasePresence.contains(CompositionalPhase::Oil);
        const bool gas = currentState.phasePresence.contains(CompositionalPhase::Gas);
        const bool water = currentState.phasePresence.contains(CompositionalPhase::Water);

        const auto phaseComposition = [&](int phase) -> const auto &
        {
            if (phase == Indices::Phase::liquid)
                return currentState.liquidMoleFraction;
            if (phase == Indices::Phase::vapor)
                return currentState.vaporMoleFraction;
            return currentState.aqueousMoleFraction;
        };

        const auto writeFugacityBlock = [&](
            const auto &equationIndices,
            int referencePhase,
            int otherPhase,
            int alternateReferencePhase = -1,
            bool alternateReferenceActive = false)
        {
            const auto &referenceComposition = phaseComposition(referencePhase);
            const auto &otherComposition = phaseComposition(otherPhase);

            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                const int equation = equationIndices[c];
                const bool referenceTrace =
                    scalarValue(referenceComposition[c]) <=
                    NaturalNumerics::phaseEquilibriumTraceComposition;
                const bool otherTrace =
                    scalarValue(otherComposition[c]) <=
                    NaturalNumerics::phaseEquilibriumTraceComposition;

                // 数值：PTz flash 可能把某组分压到机器精度以下。N-1 组成参数化
                // 无法稳定表示这种边界组分；若继续强制逸度相等会产生非物理的大残差。
                // 因此 trace 组分改用活跃边界 x_i=0。
                if (otherTrace)
                {
                    residual.value[static_cast<std::size_t>(equation)] =
                        otherComposition[c];
                    continue;
                }

                if (referenceTrace)
                {
                    // 约束：若该组分在油相为 trace、但气/水相均存在，O-G 方程块已经
                    // 提供油相边界条件；O-W 方程块改为约束剩余的 G-W 化学平衡，
                    // 避免重复写入同一油相边界方程。
                    if (alternateReferenceActive && alternateReferencePhase >= 0)
                    {
                        const auto &alternateComposition =
                            phaseComposition(alternateReferencePhase);
                        if (scalarValue(alternateComposition[c]) >
                            NaturalNumerics::phaseEquilibriumTraceComposition)
                        {
                            residual.value[static_cast<std::size_t>(equation)] =
                                (currentProperties.fugacity[static_cast<std::size_t>(alternateReferencePhase)][c] -
                                 currentProperties.fugacity[static_cast<std::size_t>(otherPhase)][c]) /
                                units::bar / fugacityScalingFactor;
                            continue;
                        }
                    }

                    residual.value[static_cast<std::size_t>(equation)] =
                        referenceComposition[c];
                    continue;
                }

                residual.value[static_cast<std::size_t>(equation)] =
                    (currentProperties.fugacity[static_cast<std::size_t>(referencePhase)][c] -
                     currentProperties.fugacity[static_cast<std::size_t>(otherPhase)][c]) /
                    units::bar / fugacityScalingFactor;
            }
        };

        const auto writeInactiveBlock = [&](
            const auto &equationIndices,
            const auto &composition,
            const Scalar &saturation)
        {
            // 状态：相不存在时，N 个逸度方程没有物理意义。这里用 N-1 个组成恒等行
            // 固定消失相组成，并用最后一行施加 S_phase=0，使局部系统继续保持方阵。
            for (int component = 0;
                 component < Indices::numIndependentCompositionsPerPhase;
                 ++component)
            {
                const Scalar &value = composition[static_cast<std::size_t>(component)];
                residual.value[static_cast<std::size_t>(
                    equationIndices[static_cast<std::size_t>(component)])] =
                    value - scalarValue(value);
            }
            residual.value[static_cast<std::size_t>(equationIndices.back())] =
                saturation;
        };

        if (oil)
        {
            if (gas)
                writeFugacityBlock(Indices::Equation::fugacity,
                                   Indices::Phase::liquid, Indices::Phase::vapor);
            else
                writeInactiveBlock(Indices::Equation::fugacity,
                                   currentState.vaporMoleFraction,
                                   currentState.vaporSaturation);

            if (water)
                writeFugacityBlock(Indices::Equation::waterFugacity,
                                   Indices::Phase::liquid, Indices::Phase::water,
                                   Indices::Phase::vapor, gas);
            else
                writeInactiveBlock(Indices::Equation::waterFugacity,
                                   currentState.aqueousMoleFraction,
                                   currentState.waterSaturation);
        }
        else
        {
            // 状态：油相缺失时，第一个方程块负责移除参考油相未知量；若气/水两相
            // 仍存在，第二个方程块直接施加 G-W 逸度平衡。
            writeInactiveBlock(Indices::Equation::fugacity,
                               currentState.liquidMoleFraction,
                               currentState.liquidSaturation);

            if (gas && water)
                writeFugacityBlock(Indices::Equation::waterFugacity,
                                   Indices::Phase::vapor, Indices::Phase::water);
            else if (gas)
                writeInactiveBlock(Indices::Equation::waterFugacity,
                                   currentState.aqueousMoleFraction,
                                   currentState.waterSaturation);
            else
                writeInactiveBlock(Indices::Equation::waterFugacity,
                                   currentState.vaporMoleFraction,
                                   currentState.vaporSaturation);
        }
    }
    else
    {
    if (currentState.hydrocarbonPhaseState ==
        HydrocarbonPhaseState::TwoPhase)
    {
        const auto fugacityResidual =
            fugacityEquilibriumResidual<Indices>(
                currentProperties,
                units::bar,
                fugacityScalingFactor);

        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            residual.value[
                static_cast<std::size_t>(
                    Indices::Equation::fugacity[
                        static_cast<std::size_t>(component)])] =
                fugacityResidual[
                    static_cast<std::size_t>(component)];
        }
    }
    else if (currentState.hydrocarbonPhaseState ==
             HydrocarbonPhaseState::LiquidOnly)
    {
        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            const Scalar &inactiveVaporComposition =
                currentState.vaporMoleFraction[
                    static_cast<std::size_t>(component)];

            residual.value[
                static_cast<std::size_t>(
                    Indices::Equation::fugacity[
                        static_cast<std::size_t>(component)])] =
                inactiveVaporComposition -
                scalarValue(inactiveVaporComposition);
        }

        residual.value[
            static_cast<std::size_t>(
                Indices::Equation::fugacity.back())] =
            currentState.vaporSaturation;
    }
    else
    {
        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            // 气单相退化行约束 vapor-composition 变量，保持方程数与未知数一致。
            const Scalar &constrainedComposition =
                currentState.vaporMoleFraction[
                    static_cast<std::size_t>(component)];

            residual.value[
                static_cast<std::size_t>(
                    Indices::Equation::fugacity[
                        static_cast<std::size_t>(component)])] =
                constrainedComposition -
                scalarValue(constrainedComposition);
        }

        residual.value[
            static_cast<std::size_t>(
                Indices::Equation::fugacity.back())] =
            currentState.liquidSaturation;
    }

    }

    residual.value[
        static_cast<std::size_t>(
            Indices::Equation::volumeClosure)] =
        saturationClosureResidual<Indices>(
            currentState);

    if constexpr (Indices::hasAqueousCO2Dissolution)
    {
        residual.value[
            static_cast<std::size_t>(
                Indices::Equation::aqueousCO2Equilibrium)] =
            aqueousCO2EquilibriumResidual<Indices>(
                currentProperties,
                currentState.hydrocarbonPhaseState,
                dissolvedCO2Component,
                units::bar,
                fugacityScalingFactor);
    }

    // 非完井单元的 well-control 行保持零残差，并通过 AD 形成单位对角。
    if constexpr (Indices::hasWellUnknown)
    {
        residual.value[
            static_cast<std::size_t>(
                Indices::Equation::wellControl)] =
            currentState.wellPressure -
            scalarValue(currentState.wellPressure);
    }

    return residual;
}

/**
 * @brief 向已组装的单元 local/base residual 增加面通量与井源项。
 *
 * 只有守恒行依赖跨单元 transport 和井源；相平衡、closure 与 dummy well
 * 行保持 local/base block 原值。该拆分允许 SNES residual 与同一状态的
 * Jacobian 复用蓄积和局部方程构造，而无需缓存逐面 AD 大对象。
 */
template <class Indices, class Scalar, class FluxScalar>
void addCellConservationTransportAndWell(
    CellResidualBlock<Indices, Scalar> &residual,
    const FaceMassFlux<Indices, FluxScalar> &outwardFaceFlux,
    const std::array<Scalar, Indices::numComponents> &wellComponentMassSource,
    const Scalar &wellWaterMassSource)
{
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        const int equation = Indices::Equation::massConservation[c];
        residual.value[static_cast<std::size_t>(equation)] +=
            outwardFaceFlux.component[c] - wellComponentMassSource[c];
    }

    if constexpr (Indices::hasIndependentWaterConservation)
    {
        residual.value[static_cast<std::size_t>(
            Indices::Equation::waterConservation)] +=
            outwardFaceFlux.water - wellWaterMassSource;
    }
}

/**
 * @brief 兼容原统一入口：local/base block + transport/source 增量。
 */
template <class Indices, class Scalar>
[[nodiscard]] CellResidualBlock<Indices, Scalar>
assembleCellResidual(
    const CellState<Indices, Scalar> &currentState,
    const CellProperties<Indices, Scalar> &currentProperties,
    const CellProperties<Indices, double> &previousProperties,
    const FaceMassFlux<Indices, Scalar> &outwardFaceFlux,
    const std::array<Scalar, Indices::numComponents> &wellComponentMassSource,
    const Scalar &wellWaterMassSource,
    double cellVolume,
    double timeStep,
    double fugacityScalingFactor,
    int dissolvedCO2Component = -1,
    double rockDensity = 2650.0,
    const std::array<double, Indices::numComponents> *standardGasDensity = nullptr)
{
    auto residual = assembleCellLocalResidual<Indices, Scalar>(
        currentState,
        currentProperties,
        previousProperties,
        cellVolume,
        timeStep,
        fugacityScalingFactor,
        dissolvedCO2Component,
        rockDensity,
        standardGasDensity);

    addCellConservationTransportAndWell<Indices>(
        residual,
        outwardFaceFlux,
        wellComponentMassSource,
        wellWaterMassSource);

    return residual;
}

} // namespace MPMC
