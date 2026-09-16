/**
 * @file ad_common_integration_test.cpp
 * @brief 单元测试：验证 `ad_common_integration` 的核心语义、边界条件和回归行为。
 */
#include <ad/Evaluation.hpp>
#include <ad/Math.hpp>
#include <common/piecewise_linear_interpolation.hpp>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    using Eval = MPMC::DenseAd::DynamicEvaluation<double, 4>;
    const Eval x = Eval::createVariable(2, 0.25, 0);
    std::vector<double> xv{0.0, 1.0};
    std::vector<double> yv{0.0, 2.0};
    const Eval y = MPMC::PiecewiseLinearInterpolation::evaluate(xv, yv, x);
    if (std::abs(y.value() - 0.5) > 1e-12 || y.size() != 2 || std::abs(y.derivative(0) - 2.0) > 1e-12 || std::abs(y.derivative(1)) > 1e-12)
        return 2;

    const Eval outside = Eval::createVariable(2, -1.0, 0);
    const Eval clamped = MPMC::PiecewiseLinearInterpolation::evaluate(xv, yv, outside);
    if (clamped.size() != 2 || std::abs(clamped.value()) > 1e-12 || std::abs(clamped.derivative(0)) > 1e-12 || std::abs(clamped.derivative(1)) > 1e-12)
        return 3;

    std::cout << "Dynamic AD/common integration: PASS\n";
    return 0;
}
