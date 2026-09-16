/**
 * @file ad_test.cpp
 * @brief 单元测试：验证 `ad` 的核心语义、边界条件和回归行为。
 */
#include <ad/DynamicEvaluation.hpp>
#include <ad/Evaluation.hpp>
#include <ad/Math.hpp>
#include <test/test_utils.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{

using Fixed = MPMC::DenseAd::Evaluation<double, 3>;
using Dynamic = MPMC::DenseAd::DynamicEvaluation<double>;

template <class Eval>
Eval function(const Eval &x, const Eval &y, const Eval &z)
{
    using namespace MPMC::DenseAd;
    return x * y + sin(z) + exp(x) / sqrt(y + 4.0) + log(z + 2.0) +
           atan2(y, x + 3.0) + pow(x + 1.0, 2.0);
}

void testFixedAndDynamicAgreement()
{
    const auto fx = Fixed::createVariable(1.2, 0);
    const auto fy = Fixed::createVariable(0.8, 1);
    const auto fz = Fixed::createVariable(0.4, 2);
    const auto fixed = function(fx, fy, fz);

    const auto dx = Dynamic::createVariable(3, 1.2, 0);
    const auto dy = Dynamic::createVariable(3, 0.8, 1);
    const auto dz = Dynamic::createVariable(3, 0.4, 2);
    const auto dynamic = function(dx, dy, dz);

    MPMC::Test::requireNear(dynamic.value(), fixed.value(), 1.0e-13,
                            "Fixed/Dynamic AD value mismatch");
    for (int i = 0; i < 3; ++i)
    {
        MPMC::Test::requireNear(dynamic.derivative(i), fixed.derivative(i), 1.0e-13,
                                "Fixed/Dynamic AD derivative mismatch");
    }
}

void testAnalyticDerivative()
{
    using Eval = MPMC::DenseAd::Evaluation<double, 1>;
    using namespace MPMC::DenseAd;

    constexpr double xValue = 1.1;
    const auto x = Eval::createVariable(xValue, 0);
    const auto g = x * x * x + sin(x);

    MPMC::Test::requireNear(g.value(),
                            xValue * xValue * xValue + std::sin(xValue),
                            1.0e-13,
                            "AD function value mismatch");
    MPMC::Test::requireNear(g.derivative(0),
                            3.0 * xValue * xValue + std::cos(xValue),
                            1.0e-13,
                            "AD first derivative mismatch");
}

void testSecondDerivative()
{
    using Inner = MPMC::DenseAd::Evaluation<double, 1>;
    using Outer = MPMC::DenseAd::Evaluation<Inner, 1>;
    using namespace MPMC::DenseAd;

    constexpr double xValue = 1.1;
    const auto innerX = Inner::createVariable(xValue, 0);
    const auto x = Outer::createVariable(innerX, 0);
    const auto g = x * x * x + sin(x);

    MPMC::Test::requireNear(g.derivative(0).derivative(0),
                            6.0 * xValue - std::sin(xValue),
                            1.0e-12,
                            "AD second derivative mismatch");
}

void testDynamicDimensionGuard()
{
    const auto lhs = Dynamic::createVariable(3, 1.0, 0);
    const auto rhs = Dynamic::createVariable(2, 2.0, 0);

    bool threw = false;
    try
    {
        const auto invalid = lhs + rhs;
        (void)invalid;
    }
    catch (const std::invalid_argument &)
    {
        threw = true;
    }

    MPMC::Test::require(threw, "Dynamic AD dimension mismatch was not rejected");
}

void testPowAtZero()
{
    using Eval = MPMC::DenseAd::Evaluation<double, 1>;
    using namespace MPMC::DenseAd;

    const auto x = Eval::createVariable(0.0, 0);
    const auto y = pow(x, Eval::createConstant(2.0));

    MPMC::Test::requireNear(y.value(), 0.0, 0.0, "pow(0,2) value mismatch");
    MPMC::Test::requireNear(y.derivative(0), 0.0, 0.0,
                            "pow(0,2) derivative mismatch");
}

} // namespace

int main()
{
    testFixedAndDynamicAgreement();
    testAnalyticDerivative();
    testSecondDerivative();
    testDynamicDimensionGuard();
    testPowAtZero();

    std::cout << "AD validation: ALL PASS\n";
    return 0;
}
