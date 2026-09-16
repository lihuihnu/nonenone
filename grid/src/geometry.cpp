/**
 * @file geometry.cpp
 * @brief CpGrid 几何量计算与向量辅助的实现。
 */
#include <cpgrid/geometry.hpp>

#include <cpgrid/node.hpp>

#include <cmath>
#include <stdexcept>

namespace MPMC
{

Point cross(
    const Point &a,
    const Point &b) noexcept
{
    return Point(
        a.y() * b.z() - a.z() * b.y(),
        a.z() * b.x() - a.x() * b.z(),
        a.x() * b.y() - a.y() * b.x());
}

double triangleArea(
    const Point &a,
    const Point &b,
    const Point &c) noexcept
{
    return 0.5 *
           cross(
               b - a,
               c - a)
               .norm();
}

double tetrahedronVolume(
    const Point &a,
    const Point &b,
    const Point &c,
    const Point &d) noexcept
{
    const Point ab = b - a;
    const Point ac = c - a;
    const Point ad = d - a;

    return std::abs(
               ab.dot(cross(ac, ad))) /
           6.0;
}

Point averagePoint(
    const std::vector<Node *> &nodes)
{
    if (nodes.empty())
    {
        throw std::invalid_argument(
            "averagePoint requires at least one node.");
    }

    Point sum;

    for (const Node *node : nodes)
    {
        if (node == nullptr)
        {
            throw std::logic_error(
                "Polygon contains a null node pointer.");
        }

        sum += static_cast<const Point &>(*node);
    }

    return sum /
           static_cast<double>(nodes.size());
}

PolygonGeometry polygonGeometry(
    const std::vector<Node *> &nodes,
    const Point &referencePoint)
{
    if (nodes.size() < 3)
    {
        throw std::invalid_argument(
            "polygonGeometry requires at least three nodes.");
    }

    PolygonGeometry result;
    Point weighted;
    Point areaVector;

    for (std::size_t i = 0;
         i < nodes.size();
         ++i)
    {
        const Node *first = nodes[i];
        const Node *second =
            nodes[(i + 1) % nodes.size()];

        if (first == nullptr ||
            second == nullptr)
        {
            throw std::logic_error(
                "Polygon contains a null node pointer.");
        }

        const Point &p1 =
            static_cast<const Point &>(*first);
        const Point &p2 =
            static_cast<const Point &>(*second);

        // 保持 polygonArea()/polygonCentroid() 的三角面积计算路径。
        const double area =
            triangleArea(
                referencePoint,
                p1,
                p2);

        result.area += area;
        weighted +=
            (referencePoint + p1 + p2) *
            (area / 3.0);

        // 保持 polygonUnitNormal() 的面积向量累积路径。
        const Point a = p1 - referencePoint;
        const Point b = p2 - referencePoint;
        areaVector += cross(a, b);
    }

    if (!std::isfinite(result.area) ||
        result.area <= 0.0)
    {
        throw std::runtime_error(
            "polygonCentroid encountered zero or invalid area.");
    }

    result.centroid =
        weighted / result.area;

    const double normalLength =
        areaVector.norm();

    if (!std::isfinite(normalLength) ||
        normalLength <= 0.0)
    {
        throw std::runtime_error(
            "polygonUnitNormal encountered a degenerate polygon.");
    }

    result.unitNormal =
        areaVector / normalLength;

    return result;
}

double polygonArea(
    const std::vector<Node *> &nodes,
    const Point &referencePoint)
{
    if (nodes.size() < 3)
    {
        throw std::invalid_argument(
            "polygonArea requires at least three nodes.");
    }

    double area = 0.0;

    for (std::size_t i = 0;
         i < nodes.size();
         ++i)
    {
        const Node *first = nodes[i];
        const Node *second =
            nodes[(i + 1) % nodes.size()];

        if (first == nullptr ||
            second == nullptr)
        {
            throw std::logic_error(
                "Polygon contains a null node pointer.");
        }

        area += triangleArea(
            referencePoint,
            static_cast<const Point &>(*first),
            static_cast<const Point &>(*second));
    }

    return area;
}

Point polygonCentroid(
    const std::vector<Node *> &nodes,
    const Point &referencePoint)
{
    if (nodes.size() < 3)
    {
        throw std::invalid_argument(
            "polygonCentroid requires at least three nodes.");
    }

    Point weighted;
    double totalArea = 0.0;

    for (std::size_t i = 0;
         i < nodes.size();
         ++i)
    {
        const Node *first = nodes[i];
        const Node *second =
            nodes[(i + 1) % nodes.size()];

        if (first == nullptr ||
            second == nullptr)
        {
            throw std::logic_error(
                "Polygon contains a null node pointer.");
        }

        const Point &p1 =
            static_cast<const Point &>(*first);
        const Point &p2 =
            static_cast<const Point &>(*second);

        const double area =
            triangleArea(
                referencePoint,
                p1,
                p2);

        weighted +=
            (referencePoint + p1 + p2) *
            (area / 3.0);

        totalArea += area;
    }

    if (!std::isfinite(totalArea) ||
        totalArea <= 0.0)
    {
        throw std::runtime_error(
            "polygonCentroid encountered zero or invalid area.");
    }

    return weighted / totalArea;
}

Point polygonUnitNormal(
    const std::vector<Node *> &nodes,
    const Point &referencePoint)
{
    if (nodes.size() < 3)
    {
        throw std::invalid_argument(
            "polygonUnitNormal requires at least three nodes.");
    }

    Point areaVector;

    for (std::size_t i = 0;
         i < nodes.size();
         ++i)
    {
        const Node *first = nodes[i];
        const Node *second =
            nodes[(i + 1) % nodes.size()];

        if (first == nullptr ||
            second == nullptr)
        {
            throw std::logic_error(
                "Polygon contains a null node pointer.");
        }

        const Point a =
            static_cast<const Point &>(*first) -
            referencePoint;
        const Point b =
            static_cast<const Point &>(*second) -
            referencePoint;

        areaVector += cross(a, b);
    }

    const double length = areaVector.norm();

    if (!std::isfinite(length) ||
        length <= 0.0)
    {
        throw std::runtime_error(
            "polygonUnitNormal encountered a degenerate polygon.");
    }

    return areaVector / length;
}

} // namespace MPMC
