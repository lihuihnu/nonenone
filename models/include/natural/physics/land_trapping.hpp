/**
 * @file land_trapping.hpp
 * @brief Land 滞留气模型及扫描过程历史状态更新。
 */
#pragma once

#include <ad/Math.hpp>
#include <common/math.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief Land 气体滞留模型的代数形式。
 *
 * 对历史最大气相饱和度 `S_g,max`，Land 关系给出最终残余气饱和度
 * `S_gr = S_g,max / (1 + C_L S_g,max)`，其中 `C_L>0` 为 Land 常数。
 * 在回吸过程中，当前总气饱和度被拆分为自由气和滞留气：
 * `S_g = S_g,free + S_g,trap`。下面的实现对 Land 曲线进行代数反演，
 * 由当前 `S_g` 和历史 `S_g,max` 直接得到 `S_g,trap`，因此无需额外
 * 非线性未知量。排驱阶段 `S_g >= S_g,max` 时滞留量为零。
 */
class LandTrappingModel final
{
public:
    explicit LandTrappingModel(double landConstant = 2.0)
        : landConstant_(landConstant)
    {
        if (!(landConstant_ > 0.0) || !std::isfinite(landConstant_))
            throw std::invalid_argument("Land constant must be finite and positive.");
    }

    /** @brief 由历史最大气饱和度计算 Land 最终残余气饱和度。 */
    [[nodiscard]] double residualGasAtMaximum(double maximumGasSaturation) const
    {
        const double sgi = std::max(maximumGasSaturation, 1.0e-12);
        return sgi / (1.0 + landConstant_ * sgi);
    }

    /**
     * @brief 返回当前时刻不可动的 Land 滞留气饱和度。
     *
     * AD 标量可直接传入，因此该代数关系会自动贡献 Newton Jacobian。
     */
    template <class Scalar>
    [[nodiscard]] Scalar trappedGasSaturation(
        const Scalar &totalGasSaturation,
        double maximumGasSaturation) const
    {
        if (scalarValue(totalGasSaturation) + 1.0e-12 >= maximumGasSaturation)
            return Scalar(0.0);

        const double sgr = residualGasAtMaximum(maximumGasSaturation);
        Scalar delta = max(totalGasSaturation - Scalar(sgr), Scalar(0.0));
        Scalar discriminant = max(
            delta * delta + (4.0 / landConstant_) * delta,
            Scalar(0.0));
        Scalar freeGas = 0.5 * (delta + sqrt(discriminant));
        Scalar nonnegativeGas = max(totalGasSaturation, Scalar(0.0));
        freeGas = min(max(freeGas, Scalar(0.0)), nonnegativeGas);
        return nonnegativeGas - freeGas;
    }

private:
    double landConstant_{2.0};
};

} // namespace MPMC
