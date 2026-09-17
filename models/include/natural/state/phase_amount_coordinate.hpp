/**
 * @file phase_amount_coordinate.hpp
 * @brief 连续的油相 phase-component amount 坐标辅助函数。
 */
#pragma once

#include <common/math.hpp>
#include <natural/numerics.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace MPMC
{

/**
 * @brief 补全 `q_i = S_o x_i` 的 dependent 末组分 `q_N`。
 *
 * 与 `x_N = 1-sum(x_i)` 相比，phase amount 坐标在 `S_o -> 0` 时自然收敛到
 * 零向量。这里保留 AD 导数，仅把机器精度内的 dependent amount 标量值平移
 * 到零，避免 `S_o-sum(q_i)` 的浮点消去伪影。
 */
template <class Scalar, std::size_t N>
void completeOilPhaseAmount(
    const Scalar &oilSaturation,
    std::array<Scalar, N> &amount)
{
    static_assert(N > 0, "Oil phase amount requires at least one component.");

    Scalar dependent = oilSaturation;
    for (std::size_t component = 0; component + 1 < N; ++component)
        dependent -= amount[component];

    const double dependentValue = scalarValue(dependent);
    const double cancellationBoundary =
        NaturalNumerics::dependentCompositionCancellationBoundary *
        std::max(1.0, std::abs(scalarValue(oilSaturation)));
    if (std::abs(dependentValue) <= cancellationBoundary)
        dependent += -dependentValue;

    amount.back() = dependent;
}

/**
 * @brief 从 `q_i=S_o x_i` 恢复活动油相的摩尔分数。
 *
 * Oil 不活动或饱和度已经退化到机器归一化尺度以下时，`x=q/S_o` 没有物理
 * 意义；此时返回调用者提供的有限 fallback。inactive residual 使用 q 本身闭合，
 * 因而 fallback 不承担任何消失相约束。
 */
template <class Scalar, std::size_t N>
[[nodiscard]] std::array<Scalar, N> oilMoleFractionFromAmount(
    const Scalar &oilSaturation,
    const std::array<Scalar, N> &amount,
    const std::array<double, N> &fallback)
{
    std::array<Scalar, N> composition{};
    if (scalarValue(oilSaturation) >
        NaturalNumerics::minimumNormalizationDenominator)
    {
        for (std::size_t component = 0; component < N; ++component)
            composition[component] = amount[component] / oilSaturation;
        return composition;
    }

    for (std::size_t component = 0; component < N; ++component)
        composition[component] = Scalar(fallback[component]);
    return composition;
}

} // namespace MPMC
