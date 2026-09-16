/**
 * @file competitive_adsorption.hpp
 * @brief 多组分竞争 Langmuir 吸附模型。
 */
#pragma once

#include <common/math.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 多组分竞争 Langmuir 吸附模型。
 *
 * theta_i = thetaMax_i b_i p_i / (1 + sum_j b_j p_j), p_i=y_i p_g.
 */
template <int NumComponents>
class CompetitiveLangmuirAdsorption final
{
public:
    using ParameterArray = std::array<double, NumComponents>;

    CompetitiveLangmuirAdsorption() = default;

    CompetitiveLangmuirAdsorption(ParameterArray thetaMax, ParameterArray coefficient)
        : thetaMax_(thetaMax), coefficient_(coefficient)
    {
        for (int i = 0; i < NumComponents; ++i)
        {
            if (!std::isfinite(thetaMax_[i]) || thetaMax_[i] < 0.0 ||
                !std::isfinite(coefficient_[i]) || coefficient_[i] < 0.0)
                throw std::invalid_argument("Langmuir parameters must be finite and nonnegative.");
        }
    }

    template <class Scalar>
    [[nodiscard]] std::array<Scalar, NumComponents>
    adsorbedVolume(
        const Scalar &gasPressure,
        const std::array<Scalar, NumComponents> &gasMoleFractions) const
    {
        std::array<Scalar, NumComponents> partialPressure{};
        Scalar denominator = 1.0;
        for (int i = 0; i < NumComponents; ++i)
        {
            partialPressure[i] = gasMoleFractions[i] * gasPressure;
            denominator += coefficient_[i] * partialPressure[i];
        }
        std::array<Scalar, NumComponents> theta{};
        for (int i = 0; i < NumComponents; ++i)
            theta[i] = thetaMax_[i] * coefficient_[i] * partialPressure[i] / denominator;
        return theta;
    }

private:
    ParameterArray thetaMax_{};
    ParameterArray coefficient_{};
};

} // namespace MPMC
