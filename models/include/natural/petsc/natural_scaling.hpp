/**
 * @file natural_scaling.hpp
 * @brief Natural 非线性系统的一致残差/变量尺度化参数。
 */
#pragma once

#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 物理变量/方程的对角尺度。
 *
 * 变量侧采用 `U = D_u \hat U`；方程侧采用 `\hat R = D_r R`。
 * 因而线性系统严格为 `D_r J D_u d\hat U = -D_r R`，随后 Newton step
 * 再乘回 `D_u`。该变换只改变数值尺度，不改变原非线性零点。
 */
struct NaturalScalingOptions final
{
    bool enabled{false};
    double pressureScale{1.0e7};       ///< Pa，对应约 100 bar。
    double compositionScale{1.0};
    double saturationScale{1.0};
    double massResidualScale{1.0};     ///< kg/s characteristic scale。
    double fugacityResidualScale{1.0}; ///< 生产 residual 已按 bar 无量纲化。
    double closureResidualScale{1.0};
    double rateWellResidualFloor{1.0}; ///< m3/s，用于 RATE 行的 characteristic scale 下限。

    void validate() const
    {
        const auto positiveFinite = [](double value) {
            return value > 0.0 && std::isfinite(value);
        };
        if (!positiveFinite(pressureScale) ||
            !positiveFinite(compositionScale) ||
            !positiveFinite(saturationScale) ||
            !positiveFinite(massResidualScale) ||
            !positiveFinite(fugacityResidualScale) ||
            !positiveFinite(closureResidualScale) ||
            !positiveFinite(rateWellResidualFloor))
        {
            throw std::invalid_argument("Natural scaling factors must be finite and positive.");
        }
    }
};

} // namespace MPMC
