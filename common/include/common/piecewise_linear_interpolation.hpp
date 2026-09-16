/**
 * @file piecewise_linear_interpolation.hpp
 * @brief 一维分段线性插值工具。
 */
#pragma once

#include <common/interval_search.hpp>

#include <cstddef>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 一维分段线性插值。
 *
 * 支持严格升序或严格降序坐标，也支持普通浮点和 AD 查询值。
 */
class PiecewiseLinearInterpolation final
{
  public:
    /**
     * @brief 计算分段线性插值值。
     *
     * 默认边界行为为 Clamp。
     */
    template <class Values, class Evaluation>
    [[nodiscard]] static Evaluation evaluate(
        const Values &xValues,
        const Values &yValues,
        const Evaluation &x,
        BoundaryPolicy boundaryPolicy = BoundaryPolicy::Clamp)
    {
        validate_(xValues, yValues);

        if (boundaryPolicy == BoundaryPolicy::Clamp)
        {
            if (xValues.front() < xValues.back())
            {
                if (scalarValue(x) <= xValues.front())
                {
                    return x * 0.0 + yValues.front();
                }

                if (scalarValue(x) >= xValues.back())
                {
                    return x * 0.0 + yValues.back();
                }
            }
            else
            {
                if (scalarValue(x) >= xValues.front())
                {
                    return x * 0.0 + yValues.front();
                }

                if (scalarValue(x) <= xValues.back())
                {
                    return x * 0.0 + yValues.back();
                }
            }
        }

        return evaluateSegment(
            xValues,
            yValues,
            x,
            findInterval(xValues, x));
    }

    /**
     * @brief 在指定线段上执行线性插值/外推。
     */
    template <class Values, class Evaluation>
    [[nodiscard]] static Evaluation evaluateSegment(
        const Values &xValues,
        const Values &yValues,
        const Evaluation &x,
        std::size_t segment)
    {
        validate_(xValues, yValues);

        if (segment + 1 >= xValues.size())
        {
            throw std::out_of_range(
                "Interpolation segment is outside the table.");
        }

        const auto x0 = xValues[segment];
        const auto x1 = xValues[segment + 1];
        const auto y0 = yValues[segment];
        const auto y1 = yValues[segment + 1];

        if (x1 == x0)
        {
            throw std::invalid_argument(
                "Interpolation segment has duplicate x coordinates.");
        }

        const auto slope =
            (y1 - y0) / (x1 - x0);

        return Evaluation(
            y0 + (x - x0) * slope);
    }

    /**
     * @brief 返回当前分段的斜率。
     *
     * Clamp 模式在表格范围外返回 0；线性外推模式返回首/末线段斜率。
     */
    template <class Values, class Evaluation>
    [[nodiscard]] static auto derivative(
        const Values &xValues,
        const Values &yValues,
        const Evaluation &x,
        BoundaryPolicy boundaryPolicy = BoundaryPolicy::Clamp)
    {
        validate_(xValues, yValues);

        if (boundaryPolicy == BoundaryPolicy::Clamp)
        {
            const auto scalar =
                scalarValue(x);

            if (xValues.front() < xValues.back())
            {
                if (scalar <= xValues.front() ||
                    scalar >= xValues.back())
                {
                    return (yValues.front() - yValues.front());
                }
            }
            else
            {
                if (scalar >= xValues.front() ||
                    scalar <= xValues.back())
                {
                    return (yValues.front() - yValues.front());
                }
            }
        }

        const std::size_t segment =
            findInterval(xValues, x);

        const auto dx =
            xValues[segment + 1] -
            xValues[segment];

        if (dx == 0)
        {
            throw std::invalid_argument(
                "Interpolation segment has duplicate x coordinates.");
        }

        return (
            yValues[segment + 1] -
            yValues[segment]) / dx;
    }

  private:
    template <class Values>
    static void validate_(
        const Values &xValues,
        const Values &yValues)
    {
        if (xValues.size() != yValues.size())
        {
            throw std::invalid_argument(
                "x/y interpolation arrays must have the same size.");
        }

        if (xValues.size() < 2)
        {
            throw std::invalid_argument(
                "Piecewise interpolation requires at least two points.");
        }

        if (!isStrictlyAscending(xValues) &&
            !isStrictlyDescending(xValues))
        {
            throw std::invalid_argument(
                "Interpolation coordinates must be strictly monotonic.");
        }
    }
};

} // namespace MPMC
