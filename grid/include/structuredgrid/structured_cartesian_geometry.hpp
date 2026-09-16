/**
 * @file structured_cartesian_geometry.hpp
 * @brief 规则笛卡尔网格可推导几何量的轻量计算工具。
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <structuredgrid/gridenums.hpp>

namespace MPMC::detail
{

/** @brief 计算单元 bulk volume。 */
[[nodiscard]] inline double structuredCellVolume(
    const std::array<double, 3> &widths) noexcept
{
    return widths[0] * widths[1] * widths[2];
}

/** @brief 计算指定面的单元中心到面中心半程距离。 */
[[nodiscard]] inline double structuredFaceHalfDistance(
    const std::array<double, 3> &widths,
    StructuredFaceOrder direction)
{
    const int axis = structuredFaceAxis(direction);
    if (axis < 0)
        throw std::invalid_argument("Invalid StructuredFace direction.");
    return 0.5 * widths[static_cast<std::size_t>(axis)];
}

/**
 * @brief 计算指定面沿其主轴的带符号面积法向分量。
 *
 * 结果与 G6 以前 `areaNormal` Vec 中对应主轴分量保持相同符号和乘法顺序。
 */
[[nodiscard]] inline double structuredFaceNormalComponent(
    const std::array<double, 3> &widths,
    StructuredFaceOrder direction)
{
    switch (direction)
    {
    case StructuredFaceOrder::XMinus: return -widths[1] * widths[2];
    case StructuredFaceOrder::XPlus:  return  widths[1] * widths[2];
    case StructuredFaceOrder::YMinus: return -widths[0] * widths[2];
    case StructuredFaceOrder::YPlus:  return  widths[0] * widths[2];
    case StructuredFaceOrder::ZMinus: return -widths[0] * widths[1];
    case StructuredFaceOrder::ZPlus:  return  widths[0] * widths[1];
    case StructuredFaceOrder::Count:
        throw std::invalid_argument("Invalid StructuredFace direction.");
    }
    throw std::invalid_argument("Invalid StructuredFace direction.");
}

/** @brief 按既有累积顺序计算一维单元中心坐标。 */
[[nodiscard]] inline std::vector<double> structuredAxisCenters(
    const std::vector<double> &widths)
{
    std::vector<double> centers(widths.size(), 0.0);
    double position = 0.0;
    for (std::size_t i = 0; i < widths.size(); ++i)
    {
        centers[i] = position + 0.5 * widths[i];
        position += widths[i];
    }
    return centers;
}

} // namespace MPMC::detail
