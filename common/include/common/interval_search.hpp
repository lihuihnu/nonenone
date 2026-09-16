/**
 * @file interval_search.hpp
 * @brief 有序区间定位与边界处理工具。
 */
#pragma once

#include <common/math.hpp>

#include <cstddef>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 一维有序表格的边界处理策略。
 */
enum class BoundaryPolicy
{
    /**
     * 超出表格范围时保持端点值。
     */
    Clamp,

    /**
     * 超出表格范围时使用首/末线段继续线性外推。
     */
    LinearExtrapolation
};

/**
 * @brief 判断坐标数组是否严格升序。
 */
template <class Values>
[[nodiscard]] bool isStrictlyAscending(const Values &values)
{
    if (values.size() < 2)
    {
        return false;
    }

    for (std::size_t i = 1; i < values.size(); ++i)
    {
        if (!(values[i - 1] < values[i]))
        {
            return false;
        }
    }

    return true;
}

/**
 * @brief 判断坐标数组是否严格降序。
 */
template <class Values>
[[nodiscard]] bool isStrictlyDescending(const Values &values)
{
    if (values.size() < 2)
    {
        return false;
    }

    for (std::size_t i = 1; i < values.size(); ++i)
    {
        if (!(values[i - 1] > values[i]))
        {
            return false;
        }
    }

    return true;
}

/**
 * @brief 在严格单调坐标表中查找查询点对应线段的左端索引。
 *
 * 无论查询值是否越界，本函数都返回合法线段索引 `[0, size-2]`。
 * 因而不会出现原降序实现返回 `size-1`、随后访问 `seg+1` 越界的问题。
 *
 * @tparam Values 坐标容器类型。
 * @tparam Evaluation 查询值类型，可以是普通浮点数或 AD 类型。
 * @param[in] values 严格升序或严格降序坐标数组。
 * @param[in] query 查询值。
 * @return 合法线段左端索引。
 * @throws std::invalid_argument 坐标数少于 2 或坐标并非严格单调。
 */
template <class Values, class Evaluation>
[[nodiscard]] std::size_t findInterval(
    const Values &values,
    const Evaluation &query)
{
    if (values.size() < 2)
    {
        throw std::invalid_argument(
            "Interpolation requires at least two coordinates.");
    }

    const bool ascending =
        values.front() < values.back();

    const bool descending =
        values.front() > values.back();

    if (!ascending && !descending)
    {
        throw std::invalid_argument(
            "Interpolation coordinates must not have identical endpoints.");
    }

    const auto x =
        scalarValue(query);

    const std::size_t last =
        values.size() - 1;

    if (ascending)
    {
        if (x <= values.front())
        {
            return 0;
        }

        if (x >= values.back())
        {
            return last - 1;
        }

        std::size_t low = 0;
        std::size_t high = last;

        while (low + 1 < high)
        {
            const std::size_t mid =
                low + (high - low) / 2;

            if (values[mid] <= x)
            {
                low = mid;
            }
            else
            {
                high = mid;
            }
        }

        return low;
    }

    if (x >= values.front())
    {
        return 0;
    }

    if (x <= values.back())
    {
        return last - 1;
    }

    std::size_t low = 0;
    std::size_t high = last;

    while (low + 1 < high)
    {
        const std::size_t mid =
            low + (high - low) / 2;

        if (values[mid] >= x)
        {
            low = mid;
        }
        else
        {
            high = mid;
        }
    }

    return low;
}

} // namespace MPMC
