/**
 * @file test_utils.hpp
 * @brief 测试源码：`test_utils`。
 */
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

namespace MPMC::Test
{

/**
 * @brief 测试条件检查。条件不满足时直接抛出异常并终止当前测试。
 */
inline void require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

/**
 * @brief 浮点数绝对误差比较。
 */
inline void requireNear(double actual,
                        double expected,
                        double tolerance,
                        const std::string &message)
{
    if (std::abs(actual - expected) > tolerance)
    {
        throw std::runtime_error(
            message + ": actual=" + std::to_string(actual) +
            ", expected=" + std::to_string(expected));
    }
}

} // namespace MPMC::Test
