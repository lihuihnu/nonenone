/**
 * @file bilinear_interpolation.hpp
 * @brief 二维双线性插值工具。
 */
#pragma once

#include <common/interval_search.hpp>

#include <cstddef>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 二维规则表格上的双线性插值。
 *
 * 表格布局固定为 `values[yIndex][xIndex]`。
 */
class BilinearInterpolation final
{
  public:
    /**
     * @brief 计算二维双线性插值。
     *
     * 默认边界行为为：
     * LinearExtrapolation。
     */
    template <
        class XValues,
        class YValues,
        class Value2D,
        class Evaluation>
    [[nodiscard]] static Evaluation evaluate(
        const XValues &xValues,
        const YValues &yValues,
        const Value2D &values,
        const Evaluation &x,
        const Evaluation &y,
        BoundaryPolicy boundaryPolicy =
            BoundaryPolicy::LinearExtrapolation)
    {
        validate_(
            xValues,
            yValues,
            values);

        const Evaluation xQuery =
            applyBoundary_(
                xValues,
                x,
                boundaryPolicy);

        const Evaluation yQuery =
            applyBoundary_(
                yValues,
                y,
                boundaryPolicy);

        const std::size_t xSegment =
            findInterval(
                xValues,
                xQuery);

        const std::size_t ySegment =
            findInterval(
                yValues,
                yQuery);

        const auto x0 =
            xValues[xSegment];
        const auto x1 =
            xValues[xSegment + 1];
        const auto y0 =
            yValues[ySegment];
        const auto y1 =
            yValues[ySegment + 1];

        if (x1 == x0 || y1 == y0)
        {
            throw std::invalid_argument(
                "Bilinear interpolation encountered a zero-width interval.");
        }

        const Evaluation tx =
            (xQuery - x0) /
            (x1 - x0);

        const Evaluation ty =
            (yQuery - y0) /
            (y1 - y0);

        const auto &f00 =
            values[ySegment][xSegment];
        const auto &f10 =
            values[ySegment][xSegment + 1];
        const auto &f01 =
            values[ySegment + 1][xSegment];
        const auto &f11 =
            values[ySegment + 1][xSegment + 1];

        return
            f00 * (1 - tx) * (1 - ty) +
            f10 * tx * (1 - ty) +
            f01 * (1 - tx) * ty +
            f11 * tx * ty;
    }

  private:
    template <class Values, class Evaluation>
    [[nodiscard]] static Evaluation applyBoundary_(
        const Values &coordinates,
        const Evaluation &query,
        BoundaryPolicy policy)
    {
        if (policy ==
            BoundaryPolicy::LinearExtrapolation)
        {
            return query;
        }

        const auto scalar =
            scalarValue(query);

        if (coordinates.front() <
            coordinates.back())
        {
            if (scalar <= coordinates.front())
            {
                return query * 0.0 +
                       coordinates.front();
            }

            if (scalar >= coordinates.back())
            {
                return query * 0.0 +
                       coordinates.back();
            }
        }
        else
        {
            if (scalar >= coordinates.front())
            {
                return query * 0.0 +
                       coordinates.front();
            }

            if (scalar <= coordinates.back())
            {
                return query * 0.0 +
                       coordinates.back();
            }
        }

        return query;
    }

    template <
        class XValues,
        class YValues,
        class Value2D>
    static void validate_(
        const XValues &xValues,
        const YValues &yValues,
        const Value2D &values)
    {
        if (xValues.size() < 2 ||
            yValues.size() < 2)
        {
            throw std::invalid_argument(
                "Bilinear interpolation requires at least 2x2 coordinates.");
        }

        if ((!isStrictlyAscending(xValues) &&
             !isStrictlyDescending(xValues)) ||
            (!isStrictlyAscending(yValues) &&
             !isStrictlyDescending(yValues)))
        {
            throw std::invalid_argument(
                "Bilinear interpolation coordinates must be strictly monotonic.");
        }

        if (values.size() != yValues.size())
        {
            throw std::invalid_argument(
                "Bilinear table row count must match y coordinates.");
        }

        for (const auto &row : values)
        {
            if (row.size() != xValues.size())
            {
                throw std::invalid_argument(
                    "Every bilinear table row must match x coordinate count.");
            }
        }
    }
};

} // namespace MPMC
