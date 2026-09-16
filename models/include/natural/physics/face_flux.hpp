/**
 * @file face_flux.hpp
 * @brief Darcy/TPFA 面通量及组分对流通量计算。
 */
#pragma once

#include <common/math.hpp>
#include <natural/numerics.hpp>
#include <natural/state/cell_state.hpp>

#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace MPMC
{

/** @brief 上游单元选择。 */
enum class UpwindSide
{
    Interior,
    Exterior
};

/**
 * @brief 单个面的质量流率结果。
 *
 * 符号约定：该值以“从 interior 单元流出为正”，
 * 面贡献，残差装配时直接累加到 interior 控制体。
 */
template <class Indices, class Scalar>
struct FaceMassFlux
{
    std::array<Scalar, Indices::numComponents> component{};
    Scalar water{0.0};
};

/**
 * @brief 一条有向连接在各相上的上游选择。
 *
 * 该对象只保存 `numPhases` 个枚举，不保存逐面 AD 值或导数。PETSc residual
 * 可用纯标量 face kernel 生成它，随后同一 CellCache generation 的 Jacobian
 * 复用完全相同的离散分支，避免 residual 为了得到数值而传播整套 AD。
 */
template <class Indices>
struct FaceFluxLinearizationDecision
{
    std::array<UpwindSide, Indices::numPhases> upwind{};
};


/**
 * @brief 重力修正后的相势差和上游判断。
 *
 * 面密度与势差按以下规则离散：
 * - 烃相面密度按两侧相饱和度加权，水相采用算术平均；
 * - 相势差 `DeltaPhi = p_ex - p_in + rho_f g (z_ex-z_in)`；
 * - 外侧压力/物性冻结 AD 导数，只对当前单元形成 Jacobian 块；
 * - 势差严格为零时按 cell volume、global id 做确定性 tie-break。
 *
 * 返回的 `UpwindSide` 决定后续 Darcy 流率中的上游流度 `lambda=k_r/mu`。
 */
[[nodiscard]] inline UpwindSide chooseUpwindSide(
    double potentialDifference,
    double interiorVolume,
    double exteriorVolume,
    long long interiorGlobalId,
    long long exteriorGlobalId) noexcept
{
    if (potentialDifference > 0.0)
        return UpwindSide::Exterior;

    if (potentialDifference < 0.0)
        return UpwindSide::Interior;

    if (interiorVolume > exteriorVolume)
        return UpwindSide::Interior;

    if (exteriorVolume > interiorVolume)
        return UpwindSide::Exterior;

    return interiorGlobalId < exteriorGlobalId
        ? UpwindSide::Interior
        : UpwindSide::Exterior;
}

template <class Indices,
          class Scalar>
[[nodiscard]] Scalar
phasePotentialDifferenceValue(
    const CellState<Indices, Scalar> &interiorState,
    const CellProperties<Indices, Scalar> &interior,
    const CellState<Indices, Scalar> &exteriorState,
    const CellProperties<Indices, Scalar> &exterior,
    int phase,
    double gravityDz)
{
    const std::size_t phaseIndex =
        static_cast<std::size_t>(phase);

    // interior 量保留 AD 引用；neighbor 只取标量值即可表达“冻结外侧导数”。
    // 避免为每个 face/phase 把外侧常量重新构造成带 N 个零导数槽的 AD 对象。
    const Scalar &rhoIn = interior.density[phaseIndex];
    const double rhoEx = scalarValue(exterior.density[phaseIndex]);

    const Scalar &satIn = interior.saturation[phaseIndex];
    const double satEx = scalarValue(exterior.saturation[phaseIndex]);

    const bool useArithmeticWaterDensity =
        Indices::hasIndependentWaterConservation &&
        phase == Indices::Phase::water;

    const Scalar averageDensity = [&]() -> Scalar
    {
        if (useArithmeticWaterDensity)
            return 0.5 * (rhoIn + rhoEx);

        const Scalar saturationAverage = 0.5 * (satIn + satEx);
        return 0.5 * (rhoIn * satIn + rhoEx * satEx) /
            std::max(scalarValue(saturationAverage), 1.0e-8);
    }();

    return scalarValue(exteriorState.pressure) +
        averageDensity * gravityDz -
        interiorState.pressure;
}

template <class Indices,
          class Scalar>
[[nodiscard]] std::pair<Scalar, UpwindSide>
phasePotentialDifference(
    const CellState<Indices, Scalar> &interiorState,
    const CellProperties<Indices, Scalar> &interior,
    const CellState<Indices, Scalar> &exteriorState,
    const CellProperties<Indices, Scalar> &exterior,
    int phase,
    double gravityDz,
    double interiorVolume,
    double exteriorVolume,
    long long interiorGlobalId,
    long long exteriorGlobalId)
{
    const Scalar difference =
        phasePotentialDifferenceValue<Indices, Scalar>(
            interiorState,
            interior,
            exteriorState,
            exterior,
            phase,
            gravityDz);

    return {
        difference,
        chooseUpwindSide(
            scalarValue(difference),
            interiorVolume,
            exteriorVolume,
            interiorGlobalId,
            exteriorGlobalId)};
}

/**
 * @brief 标准 TPFA + upstream mobility 的面总质量流率。
 *
 * 对相 `alpha`，体积流率与质量流率分别为
 * `q_alpha = -T lambda_alpha DeltaPhi_alpha`，
 * `m_dot_alpha = rho_alpha q_alpha`，其中 `lambda_alpha=k_r,alpha/mu_alpha`。
 * `T` 为两点通量近似的面传导率。
 *
 * 对水相 CO2 溶解：总水相质量流率按上游 `X_w,CO2` 拆成 CO2 与 H2O。
 */
template <class Indices,
          class Scalar>
[[nodiscard]] FaceMassFlux<Indices, Scalar>
computeFaceMassFlux(
    const CellState<Indices, Scalar> &interiorState,
    const CellProperties<Indices, Scalar> &interior,
    const CellState<Indices, Scalar> &exteriorState,
    const CellProperties<Indices, Scalar> &exterior,
    double transmissibility,
    double gravityDz,
    double interiorVolume,
    double exteriorVolume,
    long long interiorGlobalId,
    long long exteriorGlobalId,
    int dissolvedCO2Component = -1,
    const FaceFluxLinearizationDecision<Indices> *fixedDecision = nullptr)
{
    if (transmissibility < 0.0 ||
        !std::isfinite(transmissibility))
    {
        throw std::invalid_argument(
            "Transmissibility must be finite and nonnegative.");
    }

    FaceMassFlux<Indices, Scalar> result;

    for (int phase = 0;
         phase < Indices::numPhases;
         ++phase)
    {
        const Scalar potentialDifference =
            phasePotentialDifferenceValue<Indices, Scalar>(
                interiorState,
                interior,
                exteriorState,
                exterior,
                phase,
                gravityDz);

        const UpwindSide side = fixedDecision != nullptr
            ? fixedDecision->upwind[static_cast<std::size_t>(phase)]
            : chooseUpwindSide(
                scalarValue(potentialDifference),
                interiorVolume,
                exteriorVolume,
                interiorGlobalId,
                exteriorGlobalId);

        const bool interiorUpwind =
            side == UpwindSide::Interior;

        const std::size_t phaseIndex =
            static_cast<std::size_t>(phase);

        // 上游若在 interior，直接引用 AD 物性；若在 exterior，只把标量值
        // 乘入 interior-linearized 势差。两条路径都避免构造临时 frozen AD。
        const Scalar totalMassRate = [&]() -> Scalar
        {
            if (interiorUpwind)
            {
                const Scalar darcyVolumeRate =
                    -transmissibility *
                    potentialDifference *
                    interior.mobility[phaseIndex];
                return darcyVolumeRate * interior.density[phaseIndex];
            }

            const double mobility =
                scalarValue(exterior.mobility[phaseIndex]);
            const double density =
                scalarValue(exterior.density[phaseIndex]);
            const Scalar darcyVolumeRate =
                -transmissibility *
                potentialDifference *
                mobility;
            return darcyVolumeRate * density;
        }();

        if constexpr (Indices::hasIndependentWaterConservation)
        {
            if (phase == Indices::Phase::water)
            {
                if constexpr (Indices::hasAqueousCO2Dissolution)
                {
                    if (dissolvedCO2Component < 0 ||
                        dissolvedCO2Component >=
                            Indices::numComponents)
                    {
                        throw std::invalid_argument(
                            "A valid CO2 component index is required for aqueous flux.");
                    }

                    if (interiorUpwind)
                    {
                        const Scalar &aqueousCO2MassFraction =
                            interior.aqueousCO2MassFraction;
                        result.component[
                            static_cast<std::size_t>(
                                dissolvedCO2Component)] +=
                            totalMassRate * aqueousCO2MassFraction;
                        result.water +=
                            totalMassRate *
                            (1.0 - aqueousCO2MassFraction);
                    }
                    else
                    {
                        const double aqueousCO2MassFraction =
                            scalarValue(exterior.aqueousCO2MassFraction);
                        result.component[
                            static_cast<std::size_t>(
                                dissolvedCO2Component)] +=
                            totalMassRate * aqueousCO2MassFraction;
                        result.water +=
                            totalMassRate *
                            (1.0 - aqueousCO2MassFraction);
                    }
                }
                else
                {
                    result.water +=
                        totalMassRate;
                }

                continue;
            }
        }

        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            if (interiorUpwind)
            {
                result.component[c] +=
                    totalMassRate * interior.massFraction[phaseIndex][c];
            }
            else
            {
                result.component[c] +=
                    totalMassRate *
                    scalarValue(exterior.massFraction[phaseIndex][c]);
            }
        }
    }

    return result;
}

/**
 * @brief Residual 专用的纯标量 TPFA face kernel。
 *
 * 输入仍可来自 AD `CellState/CellProperties`，但这里只读取 `scalarValue()`；
 * 因而 residual 不再为每个 face/phase/component 构造和传播 N 个导数槽。
 * 可选 `decision` 记录各相上游选择，供同一状态随后 Jacobian 复用离散分支。
 */
template <class Indices, class Scalar>
[[nodiscard]] FaceMassFlux<Indices, double>
computeFaceMassFluxValue(
    const CellState<Indices, Scalar> &interiorState,
    const CellProperties<Indices, Scalar> &interior,
    const CellState<Indices, Scalar> &exteriorState,
    const CellProperties<Indices, Scalar> &exterior,
    double transmissibility,
    double gravityDz,
    double interiorVolume,
    double exteriorVolume,
    long long interiorGlobalId,
    long long exteriorGlobalId,
    int dissolvedCO2Component = -1,
    FaceFluxLinearizationDecision<Indices> *decision = nullptr)
{
    if (transmissibility < 0.0 || !std::isfinite(transmissibility))
        throw std::invalid_argument(
            "Transmissibility must be finite and nonnegative.");

    FaceMassFlux<Indices, double> result;

    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        const double rhoIn = scalarValue(interior.density[p]);
        const double rhoEx = scalarValue(exterior.density[p]);
        const double satIn = scalarValue(interior.saturation[p]);
        const double satEx = scalarValue(exterior.saturation[p]);

        const bool waterArithmetic =
            Indices::hasIndependentWaterConservation &&
            phase == Indices::Phase::water;

        const double averageDensity = waterArithmetic
            ? 0.5 * (rhoIn + rhoEx)
            : 0.5 * (rhoIn * satIn + rhoEx * satEx) /
                std::max(0.5 * (satIn + satEx), 1.0e-8);

        const double potentialDifference =
            scalarValue(exteriorState.pressure) +
            averageDensity * gravityDz -
            scalarValue(interiorState.pressure);

        const UpwindSide side = chooseUpwindSide(
            potentialDifference,
            interiorVolume,
            exteriorVolume,
            interiorGlobalId,
            exteriorGlobalId);
        if (decision != nullptr)
            decision->upwind[p] = side;

        const bool interiorUpwind = side == UpwindSide::Interior;
        const double mobility = interiorUpwind
            ? scalarValue(interior.mobility[p])
            : scalarValue(exterior.mobility[p]);
        const double density = interiorUpwind
            ? scalarValue(interior.density[p])
            : scalarValue(exterior.density[p]);
        const double totalMassRate =
            -transmissibility * potentialDifference * mobility * density;

        if constexpr (Indices::hasIndependentWaterConservation)
        {
            if (phase == Indices::Phase::water)
            {
                if constexpr (Indices::hasAqueousCO2Dissolution)
                {
                    if (dissolvedCO2Component < 0 ||
                        dissolvedCO2Component >= Indices::numComponents)
                        throw std::invalid_argument(
                            "A valid CO2 component index is required for aqueous flux.");

                    const double xco2 = interiorUpwind
                        ? scalarValue(interior.aqueousCO2MassFraction)
                        : scalarValue(exterior.aqueousCO2MassFraction);
                    result.component[static_cast<std::size_t>(dissolvedCO2Component)] +=
                        totalMassRate * xco2;
                    result.water += totalMassRate * (1.0 - xco2);
                }
                else
                {
                    result.water += totalMassRate;
                }
                continue;
            }
        }

        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            const double massFraction = interiorUpwind
                ? scalarValue(interior.massFraction[p][c])
                : scalarValue(exterior.massFraction[p][c]);
            result.component[c] += totalMassRate * massFraction;
        }
    }

    return result;
}

} // namespace MPMC
