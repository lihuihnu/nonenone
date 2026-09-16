/**
 * @file polyhedron.cpp
 * @brief CpGrid 单元多面体及输入/当前编号信息的实现。
 */
#include <cpgrid/polyhedron.hpp>

#include <cpgrid/config.hpp>
#include <cpgrid/geometry.hpp>

#include <cmath>
#include <ostream>
#include <stdexcept>

namespace MPMC
{

Polyhedron::Polyhedron(
    Polyhedron &&other) noexcept
    : currentId_(other.currentId_),
      storageIndex_(other.storageIndex_),
      faces_(std::move(other.faces_)),
      processorId_(other.processorId_),
      centroid_(other.centroid_),
      volume_(other.volume_)
{
    rebindFaceOwners();
    other.resetMovedFrom_();
}

Polyhedron &Polyhedron::operator=(
    Polyhedron &&other) noexcept
{
    if (this == &other)
    {
        return *this;
    }

    currentId_ = other.currentId_;
    storageIndex_ = other.storageIndex_;
    faces_ = std::move(other.faces_);
    processorId_ = other.processorId_;
    centroid_ = other.centroid_;
    volume_ = other.volume_;

    rebindFaceOwners();
    other.resetMovedFrom_();

    return *this;
}

void Polyhedron::resetMovedFrom_() noexcept
{
    currentId_ = 0;
    storageIndex_ = 0;
    processorId_ = 0;
    centroid_ = {};
    volume_ = 0.0;
}

void Polyhedron::rebindFaceOwners() noexcept
{
    for (Face &face : faces_)
    {
        face.setOwnerCell(this);
    }
}

Point Polyhedron::vertexAverage() const
{
    if (faces_.empty())
    {
        throw std::logic_error(
            "Polyhedron::vertexAverage requires at least one face.");
    }

    Point sum;
    std::size_t count = 0;

    /*
     * 同一个顶点会出现在多个面中。这里只需要一个稳定的内部参考点，
     * 因此允许重复计数；最终仍位于单元顶点包围区域内。
     */
    for (const Face &face : faces_)
    {
        for (const Node *node : face.nodes())
        {
            if (node == nullptr)
            {
                throw std::logic_error(
                    "Polyhedron contains a null node pointer.");
            }

            sum += static_cast<const Point &>(*node);
            ++count;
        }
    }

    if (count == 0)
    {
        throw std::logic_error(
            "Polyhedron contains no nodes.");
    }

    return sum /
           static_cast<double>(count);
}

double Polyhedron::computeVolume() const
{
#if MESH_DIM == 2
    const Point reference = vertexAverage();
    double area = 0.0;

    for (const Face &face : faces_)
    {
        if (face.nodeCount() != 2)
        {
            throw std::runtime_error(
                "A 2D Polyhedron face must contain two nodes.");
        }

        const auto &nodes = face.nodes();

        area += triangleArea(
            reference,
            static_cast<const Point &>(*nodes[0]),
            static_cast<const Point &>(*nodes[1]));
    }

    return area;
#else
    const Point cellReference = vertexAverage();
    double volume = 0.0;

    for (const Face &face : faces_)
    {
        const Point faceReference =
            face.averagePoint();

        const auto &nodes = face.nodes();

        if (nodes.size() < 3)
        {
            throw std::runtime_error(
                "A 3D Polyhedron face must contain at least three nodes.");
        }

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
                    "Polyhedron contains a null node pointer.");
            }

            volume += tetrahedronVolume(
                cellReference,
                faceReference,
                static_cast<const Point &>(*first),
                static_cast<const Point &>(*second));
        }
    }

    return volume;
#endif
}

PolyhedronGeometry Polyhedron::computeGeometry() const
{
#if MESH_DIM == 2
    const Point reference = vertexAverage();
    Point weighted;
    double totalArea = 0.0;

    for (const Face &face : faces_)
    {
        const auto &nodes = face.nodes();

        if (nodes.size() != 2 ||
            nodes[0] == nullptr ||
            nodes[1] == nullptr)
        {
            throw std::runtime_error(
                "Invalid 2D face while computing cell centroid.");
        }

        const Point &p1 =
            static_cast<const Point &>(*nodes[0]);
        const Point &p2 =
            static_cast<const Point &>(*nodes[1]);

        const double area =
            triangleArea(reference, p1, p2);

        weighted +=
            (reference + p1 + p2) *
            (area / 3.0);
        totalArea += area;
    }

    if (!std::isfinite(totalArea) ||
        totalArea <= 0.0)
    {
        throw std::runtime_error(
            "Polyhedron centroid encountered invalid area.");
    }

    return {totalArea, weighted / totalArea};
#else
    const Point cellReference = vertexAverage();
    Point weighted;
    double totalVolume = 0.0;

    for (const Face &face : faces_)
    {
        const Point faceReference =
            face.averagePoint();
        const auto &nodes = face.nodes();

        if (nodes.size() < 3)
        {
            throw std::runtime_error(
                "Invalid 3D face while computing cell centroid.");
        }

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
                    "Polyhedron contains a null node pointer.");
            }

            const Point &p1 =
                static_cast<const Point &>(*first);
            const Point &p2 =
                static_cast<const Point &>(*second);

            const double volume =
                tetrahedronVolume(
                    cellReference,
                    faceReference,
                    p1,
                    p2);

            const Point tetraCentroid =
                (cellReference +
                 faceReference +
                 p1 +
                 p2) /
                4.0;

            weighted +=
                tetraCentroid * volume;
            totalVolume += volume;
        }
    }

    if (!std::isfinite(totalVolume) ||
        totalVolume <= 0.0)
    {
        throw std::runtime_error(
            "Polyhedron centroid encountered invalid volume.");
    }

    return {totalVolume, weighted / totalVolume};
#endif
}

Point Polyhedron::centroid() const
{
#if MESH_DIM == 2
    const Point reference = vertexAverage();
    Point weighted;
    double totalArea = 0.0;

    for (const Face &face : faces_)
    {
        const auto &nodes = face.nodes();

        if (nodes.size() != 2 ||
            nodes[0] == nullptr ||
            nodes[1] == nullptr)
        {
            throw std::runtime_error(
                "Invalid 2D face while computing cell centroid.");
        }

        const Point &p1 =
            static_cast<const Point &>(*nodes[0]);
        const Point &p2 =
            static_cast<const Point &>(*nodes[1]);

        const double area =
            triangleArea(reference, p1, p2);

        weighted +=
            (reference + p1 + p2) *
            (area / 3.0);

        totalArea += area;
    }

    if (!std::isfinite(totalArea) ||
        totalArea <= 0.0)
    {
        throw std::runtime_error(
            "Polyhedron centroid encountered invalid area.");
    }

    return weighted / totalArea;
#else
    const Point cellReference = vertexAverage();
    Point weighted;
    double totalVolume = 0.0;

    for (const Face &face : faces_)
    {
        const Point faceReference =
            face.averagePoint();
        const auto &nodes = face.nodes();

        if (nodes.size() < 3)
        {
            throw std::runtime_error(
                "Invalid 3D face while computing cell centroid.");
        }

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
                    "Polyhedron contains a null node pointer.");
            }

            const Point &p1 =
                static_cast<const Point &>(*first);
            const Point &p2 =
                static_cast<const Point &>(*second);

            const double volume =
                tetrahedronVolume(
                    cellReference,
                    faceReference,
                    p1,
                    p2);

            const Point tetraCentroid =
                (cellReference +
                 faceReference +
                 p1 +
                 p2) /
                4.0;

            weighted +=
                tetraCentroid * volume;
            totalVolume += volume;
        }
    }

    if (!std::isfinite(totalVolume) ||
        totalVolume <= 0.0)
    {
        throw std::runtime_error(
            "Polyhedron centroid encountered invalid volume.");
    }

    return weighted / totalVolume;
#endif
}

std::ostream &operator<<(
    std::ostream &os,
    const Polyhedron &cell)
{
    return os
           << "Cell(id="
           << cell.currentId_
           << ", input="
           << cell.storageIndex_
           << ", ownerRank="
           << cell.processorId_
           << ", faces="
           << cell.faces_.size()
           << ')';
}

} // namespace MPMC
