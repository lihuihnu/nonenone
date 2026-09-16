/**
 * @file well_source.hpp
 * @brief 井穿孔对单元组分守恒方程的源汇项。
 */
#pragma once

#include <common/math.hpp>
#include <natural/numerics.hpp>
#include <natural/state/cell_state.hpp>
#include <well/types.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

template <class Indices, class Scalar>
struct PerforationWellResult
{
    std::array<Scalar, Indices::numPhases> surfacePhaseRate{};
    std::array<Scalar, Indices::numPhases> reservoirPhaseRate{};
    std::array<Scalar, Indices::numPhases> phaseMassRate{};
    std::array<Scalar, Indices::numComponents> componentMassSource{};
    Scalar waterMassSource{0.0};
};

/** @brief 穿孔井源项内部可复用的 AD 工作区。 */
template <class Indices, class Scalar>
struct PerforationWellWorkspace
{
    std::array<Scalar, Indices::numPhases> pressureDrop{};
    std::array<Scalar, Indices::numPhases> densityRatio{};
};

/**
 * @brief 将一个完井段的 Peaceman 型井源项直接写入调用方 scratch。
 *
 * 对相 alpha 定义 `Delta p_alpha = p_alpha - p_bhp`、
 * `b_alpha = rho_alpha/rho_alpha,sc`。采出或反向流时使用
 * `q_alpha,sc = -WI lambda_alpha Delta p_alpha b_alpha`，其中
 * `lambda_alpha = k_r,alpha/mu_alpha`。储层体积流量为
 * `q_alpha,res = q_alpha,sc/b_alpha`，质量流率为
 * `m_dot_alpha = rho_alpha q_alpha,res`。
 *
 * 注入井采用总流度 `sum(lambda_alpha)` 计算井筒总流动能力，再按
 * `injectionPhaseFraction` 分配各相；组分质量源项再按注入组分质量分数
 * 分配。统一符号约定：向油藏注入为正，从油藏采出为负。
 *
 * 与 computePerforationWellSource() 数学完全相同，但不会为每个穿孔重新
 * 默认构造整份 AD 结果对象。Natural PETSc runtime 用它复用固定 scratch；
 * 传统值返回 API 保留在下方作为兼容包装。
 */
template <class Indices, class Scalar>
void computePerforationWellSourceInto(
    PerforationWellResult<Indices, Scalar> &result,
    PerforationWellWorkspace<Indices, Scalar> &workspace,
    WellType wellType,
    double wellIndex,
    const Scalar &bottomHolePressure,
    const std::array<Scalar, Indices::numPhases> &phasePressure,
    const std::array<Scalar, Indices::numPhases> &phaseDensity,
    const std::array<Scalar, Indices::numPhases> &phaseMobility,
    const std::array<double, Indices::numPhases> &surfaceDensity,
    const std::array<double, Indices::numPhases> &injectionPhaseFraction,
    const std::array<double, Indices::numComponents> &injectionComponentMassFraction,
    const std::array<std::array<Scalar, Indices::numComponents>, Indices::numPhases> &cellMassFraction,
    const Scalar &aqueousCO2MassFraction,
    int dissolvedCO2Component)
{
    if (!(wellIndex >= 0.0) || !std::isfinite(wellIndex))
        throw std::invalid_argument("Well index must be finite and nonnegative.");

    auto &dP = workspace.pressureDrop;
    auto &b = workspace.densityRatio;
    Scalar volumeRatio = 0.0;
    Scalar totalMobility = 0.0;

    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        if (!(surfaceDensity[p] > 0.0))
            throw std::invalid_argument("Surface density must be positive.");
        dP[p] = phasePressure[p] - bottomHolePressure;
        b[p] = phaseDensity[p] / surfaceDensity[p];
        volumeRatio += injectionPhaseFraction[p] / b[p];
        totalMobility += phaseMobility[p];
    }

    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        const bool injectingThroughPhase =
            wellType == WellType::Injector && scalarValue(dP[p]) < 0.0;

        if (injectingThroughPhase)
        {
            if (!(scalarValue(volumeRatio) > 0.0))
                throw std::runtime_error("Injector volume-ratio denominator must be positive.");
            result.surfacePhaseRate[p] =
                -wellIndex * totalMobility * dP[p] / volumeRatio * injectionPhaseFraction[p];
        }
        else
        {
            result.surfacePhaseRate[p] =
                -wellIndex * phaseMobility[p] * dP[p] * b[p];
        }
        result.reservoirPhaseRate[p] = result.surfacePhaseRate[p] / b[p];
        result.phaseMassRate[p] = result.reservoirPhaseRate[p] * phaseDensity[p];
    }

    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        // 物理：三相共享同一组守恒组分。注入使用给定注入流体质量分数，
        // 采出使用井所在单元的相内质量分数；H2O 已作为普通组分守恒，
        // 因此全组分三相模式不再单独设置水质量源项。
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            Scalar source = 0.0;
            for (int phase = 0; phase < Indices::numPhases; ++phase)
            {
                const std::size_t p = static_cast<std::size_t>(phase);
                if (scalarValue(result.phaseMassRate[p]) > 0.0)
                    source += result.phaseMassRate[p] * injectionComponentMassFraction[c];
                else
                    source += result.phaseMassRate[p] * cellMassFraction[p][c];
            }
            result.componentMassSource[c] = source;
        }
        result.waterMassSource = Scalar(0.0);
        return;
    }

    if constexpr (Indices::hasIndependentWaterConservation)
    {
        const std::size_t water = static_cast<std::size_t>(Indices::Phase::water);
        const Scalar totalWaterMass = result.phaseMassRate[water];

        if constexpr (Indices::hasAqueousCO2Dissolution)
        {
            if (dissolvedCO2Component < 0 || dissolvedCO2Component >= Indices::numComponents)
                throw std::invalid_argument("Valid dissolved CO2 component index required.");

            // 注入水默认不含溶解 CO2；采出水使用单元当前质量分数。
            const Scalar sourceX = scalarValue(totalWaterMass) < 0.0
                ? aqueousCO2MassFraction : Scalar(0.0);

            result.waterMassSource = totalWaterMass * (1.0 - sourceX);
        }
        else
        {
            result.waterMassSource = totalWaterMass;
        }
    }

    else
    {
        // in-place scratch may contain a previous perforation result.
        result.waterMassSource = Scalar(0.0);
    }

    const std::size_t liquid = static_cast<std::size_t>(Indices::Phase::liquid);
    const std::size_t vapor = static_cast<std::size_t>(Indices::Phase::vapor);
    const bool injectingLiquid = scalarValue(result.phaseMassRate[liquid]) > 0.0;
    const bool injectingVapor = scalarValue(result.phaseMassRate[vapor]) > 0.0;

    const auto hydrocarbonComponentSource = [&](std::size_t component) -> Scalar
    {
        Scalar source = 0.0;
        if (injectingLiquid)
            source += result.phaseMassRate[liquid] *
                      injectionComponentMassFraction[component];
        else
            source += cellMassFraction[liquid][component] *
                      result.phaseMassRate[liquid];

        if (injectingVapor)
            source += result.phaseMassRate[vapor] *
                      injectionComponentMassFraction[component];
        else
            source += cellMassFraction[vapor][component] *
                      result.phaseMassRate[vapor];
        return source;
    };

    if constexpr (Indices::hasAqueousCO2Dissolution)
    {
        static_assert(Indices::hasWater,
                      "Aqueous CO2 dissolution requires an explicit water phase.");

        const std::size_t water = static_cast<std::size_t>(Indices::Phase::water);
        const Scalar totalWaterMass = result.phaseMassRate[water];
        const Scalar sourceX = scalarValue(totalWaterMass) < 0.0
            ? aqueousCO2MassFraction : Scalar(0.0);
        const Scalar aqueousCO2Source = totalWaterMass * sourceX;

        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            Scalar source = hydrocarbonComponentSource(c);
            if (component == dissolvedCO2Component)
                source += aqueousCO2Source;
            result.componentMassSource[c] = source;
        }
    }
    else
    {
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            result.componentMassSource[c] = hydrocarbonComponentSource(c);
        }
    }
}

/**
 * @brief 值返回兼容包装；独立调用者仍可沿用原 API。
 */
template <class Indices, class Scalar>
[[nodiscard]] PerforationWellResult<Indices, Scalar>
computePerforationWellSource(
    WellType wellType,
    double wellIndex,
    const Scalar &bottomHolePressure,
    const std::array<Scalar, Indices::numPhases> &phasePressure,
    const std::array<Scalar, Indices::numPhases> &phaseDensity,
    const std::array<Scalar, Indices::numPhases> &phaseMobility,
    const std::array<double, Indices::numPhases> &surfaceDensity,
    const std::array<double, Indices::numPhases> &injectionPhaseFraction,
    const std::array<double, Indices::numComponents> &injectionComponentMassFraction,
    const std::array<std::array<Scalar, Indices::numComponents>, Indices::numPhases> &cellMassFraction,
    [[maybe_unused]] Scalar aqueousCO2MassFraction = Scalar(0.0),
    int dissolvedCO2Component = -1)
{
    PerforationWellResult<Indices, Scalar> result;
    PerforationWellWorkspace<Indices, Scalar> workspace;
    computePerforationWellSourceInto<Indices, Scalar>(
        result,
        workspace,
        wellType,
        wellIndex,
        bottomHolePressure,
        phasePressure,
        phaseDensity,
        phaseMobility,
        surfaceDensity,
        injectionPhaseFraction,
        injectionComponentMassFraction,
        cellMassFraction,
        aqueousCO2MassFraction,
        dissolvedCO2Component);
    return result;
}

/**
 * @brief 从各相地面体积流量中提取当前 RATE/OIL/GAS/WATER 控制量。
 *
 * 返回值沿用模拟器符号：注入为正、采出为负。BHP 不是流量控制，
 * 因而传入 `WellControl::Bhp` 会抛出异常。
 */
template <class Indices, class Scalar>
[[nodiscard]] Scalar selectControlledRate(
    WellControl control,
    const std::array<Scalar, Indices::numPhases> &surfacePhaseRate,
    const std::array<Scalar, Indices::numPhases> &reservoirPhaseRate)
{
    switch (control)
    {
    case WellControl::TotalRate:
    {
        Scalar result = 0.0;
        for (const auto &rate : surfacePhaseRate)
            result += rate;
        return result;
    }
    case WellControl::ReservoirTotalRate:
    {
        Scalar result = 0.0;
        for (const auto &rate : reservoirPhaseRate)
            result += rate;
        return result;
    }
    case WellControl::OilRate:
        return surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::liquid)];
    case WellControl::GasRate:
        return surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::vapor)];
    case WellControl::WaterRate:
        if constexpr (Indices::hasWater)
            return surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::water)];
        else
            throw std::logic_error("Water-rate control requires a water phase.");
    case WellControl::Bhp:
        throw std::logic_error("BHP control is not a rate selector.");
    }
    throw std::logic_error("Unknown well control.");
}

/** Backward-compatible selector for the original surface-rate controls. */
template <class Indices, class Scalar>
[[nodiscard]] Scalar selectControlledRate(
    WellControl control,
    const std::array<Scalar, Indices::numPhases> &surfacePhaseRate)
{
    return selectControlledRate<Indices>(
        control, surfacePhaseRate, surfacePhaseRate);
}

/**
 * @brief 计算井控制方程的物理残差。
 *
 * BHP 控制：`R_w = p_bhp - p_target`；
 * 流量控制：`R_w = q_control - q_target`。
 * 这里不加入数值稳定用的 Jacobian 对角修正，该修正由装配层单独处理。
 */
template <class Scalar>
[[nodiscard]] Scalar wellControlResidual(
    WellControl control,
    const Scalar &controlledValue,
    double target,
    const Scalar &bottomHolePressure)
{
    if (control == WellControl::Bhp)
        return bottomHolePressure - target;
    return controlledValue - target;
}

} // namespace MPMC
