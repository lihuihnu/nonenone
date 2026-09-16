/**
 * @file face.cpp
 * @brief CpGrid 面拓扑、几何和相邻单元信息的实现。
 */
#include <cpgrid/face.hpp>

#include <cpgrid/config.hpp>
#include <cpgrid/geometry.hpp>
#include <cpgrid/polyhedron.hpp>

#include <cmath>
#include <ostream>
#include <stdexcept>

namespace MPMC
{

Point Face::averagePoint() const
{
    return MPMC::averagePoint(nodes_);
}

Point Face::centroid() const
{
#if MESH_DIM == 2
    if (nodes_.size() != 2 ||
        nodes_[0] == nullptr ||
        nodes_[1] == nullptr)
    {
        throw std::runtime_error(
            "A 2D face must contain exactly two valid nodes.");
    }

    return (
               static_cast<const Point &>(*nodes_[0]) +
               static_cast<const Point &>(*nodes_[1])) /
           2.0;
#else
    const Point reference =
        MPMC::averagePoint(nodes_);

    return polygonCentroid(
        nodes_,
        reference);
#endif
}

Point Face::unitNormal() const
{
#if MESH_DIM == 2
    if (nodes_.size() != 2 ||
        nodes_[0] == nullptr ||
        nodes_[1] == nullptr)
    {
        throw std::runtime_error(
            "A 2D face must contain exactly two valid nodes.");
    }

    const Point edge =
        static_cast<const Point &>(*nodes_[1]) -
        static_cast<const Point &>(*nodes_[0]);

    const double length = edge.norm();

    if (!std::isfinite(length) ||
        length <= 0.0)
    {
        throw std::runtime_error(
            "A 2D face has zero or invalid length.");
    }

    return Point(
        -edge.y() / length,
        edge.x() / length,
        0.0);
#else
    const Point reference =
        MPMC::averagePoint(nodes_);

    return polygonUnitNormal(
        nodes_,
        reference);
#endif
}

PolygonGeometry Face::computeGeometry() const
{
#if MESH_DIM == 2
    if (nodes_.size() != 2 ||
        nodes_[0] == nullptr ||
        nodes_[1] == nullptr)
    {
        throw std::runtime_error(
            "A 2D face must contain exactly two valid nodes.");
    }

    const Point &first =
        static_cast<const Point &>(*nodes_[0]);
    const Point &second =
        static_cast<const Point &>(*nodes_[1]);
    const Point edge = second - first;
    const double length = edge.norm();

    if (!std::isfinite(length) ||
        length <= 0.0)
    {
        throw std::runtime_error(
            "A 2D face has zero or invalid length.");
    }

    return {
        length,
        (first + second) / 2.0,
        Point(-edge.y() / length, edge.x() / length, 0.0)};
#else
    const Point reference =
        MPMC::averagePoint(nodes_);
    return polygonGeometry(nodes_, reference);
#endif
}

double Face::computeArea() const
{
#if MESH_DIM == 2
    if (nodes_.size() != 2 ||
        nodes_[0] == nullptr ||
        nodes_[1] == nullptr)
    {
        throw std::runtime_error(
            "A 2D face must contain exactly two valid nodes.");
    }

    return (
               static_cast<const Point &>(*nodes_[1]) -
               static_cast<const Point &>(*nodes_[0]))
        .norm();
#else
    const Point reference =
        MPMC::averagePoint(nodes_);

    return polygonArea(
        nodes_,
        reference);
#endif
}

std::ostream &operator<<(
    std::ostream &os,
    const Face &face)
{
    os << "Face(inputId="
       << face.inputFaceId_
       << ", pos="
       << face.inputPosition_
       << ", nodes="
       << face.nodes_.size()
       << ", boundary="
       << (face.isBoundary() ? "true" : "false")
       << ')';

    return os;
}

} // namespace MPMC
