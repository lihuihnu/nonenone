/**
 * @file ad_dimension_test.cpp
 * @brief 单元测试：验证 `ad_dimension` 的核心语义、边界条件和回归行为。
 */
#include <ad/Evaluation.hpp>
#include <ad/Math.hpp>
#include <common/math.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{

template <int N>
void testDimension()
{
    using Evaluation =
        MPMC::DenseAd::Evaluation<double, N>;

    const auto x =
        Evaluation::createVariable(
            N,
            2.0,
            0);

    const auto y =
        Evaluation::createVariable(
            N,
            3.0,
            N - 1);

    const Evaluation result =
        x * x +
        2.0 * x * y +
        MPMC::DenseAd::sin(y);

    const double expectedValue =
        4.0 +
        12.0 +
        std::sin(3.0);

    const double expectedDx =
        2.0 * 2.0 +
        2.0 * 3.0;

    const double expectedDy =
        2.0 * 2.0 +
        std::cos(3.0);

    const double tolerance =
        1.0e-12;

    if (std::abs(
            result.value() -
            expectedValue) >
        tolerance)
    {
        throw std::runtime_error(
            "fixed Evaluation value test failed");
    }

    if (std::abs(
            result.derivative(0) -
            expectedDx) >
        tolerance)
    {
        throw std::runtime_error(
            "fixed Evaluation first derivative test failed");
    }

    if constexpr (N > 1)
    {
        if (std::abs(
                result.derivative(N - 1) -
                expectedDy) >
            tolerance)
        {
            throw std::runtime_error(
                "fixed Evaluation last derivative test failed");
        }
    }

    for (int derivative = 1;
         derivative < N - 1;
         ++derivative)
    {
        if (std::abs(
                result.derivative(derivative)) >
            tolerance)
        {
            throw std::runtime_error(
                "fixed Evaluation inactive derivative is nonzero");
        }
    }
}

} // namespace

int main()
{
    /*
     * 覆盖原工程曾经为其维护独立 EvaluationXX.hpp 的全部固定维数。
     */
    testDimension<7>();
    testDimension<8>();
    testDimension<9>();
    testDimension<11>();
    testDimension<14>();
    testDimension<15>();
    testDimension<16>();
    testDimension<19>();
    testDimension<20>();
    testDimension<43>();

    /*
     * 同时证明新增任意维数不再需要创建 EvaluationXX.hpp。
     */
    testDimension<6>();
    testDimension<10>();
    testDimension<17>();
    testDimension<23>();

    std::cout
        << "Generic fixed-size Evaluation dimensions: ALL PASS\n";

    return 0;
}
