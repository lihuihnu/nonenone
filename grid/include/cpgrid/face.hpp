/**
 * @file face.hpp
 * @brief CpGrid 面拓扑、几何和相邻单元信息。
 */
#pragma once

#include <cpgrid/node.hpp>
#include <cpgrid/point.hpp>
#include <cpgrid/geometry.hpp>

#include <petscsys.h>

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

namespace MPMC
{

class Polyhedron;

/**
 * @brief 非结构化网格中的多边形面。
 *
 * Face 拥有自身的节点指针数组，但不拥有 Node 或 Polyhedron。
 * owner/neighbor 指针仅观察 Mesh 内稳定存储的单元。
 */
class Face final
{
  public:
    Face() = default;

    Face(
        std::vector<Node *> nodes,
        PetscInt inputPosition,
        PetscInt inputFaceId) noexcept
        : nodes_(std::move(nodes)),
          inputPosition_(inputPosition),
          inputFaceId_(inputFaceId)
    {
    }

    Face(const Face &) = default;
    Face &operator=(const Face &) = default;
    Face(Face &&) noexcept = default;
    Face &operator=(Face &&) noexcept = default;
    ~Face() = default;

    [[nodiscard]] PetscInt inputFaceId() const noexcept
    {
        return inputFaceId_;
    }

    [[nodiscard]] PetscInt inputPosition() const noexcept
    {
        return inputPosition_;
    }

    [[nodiscard]] std::size_t nodeCount() const noexcept
    {
        return nodes_.size();
    }

    [[nodiscard]] const std::vector<Node *> &
    nodes() const noexcept
    {
        return nodes_;
    }

    [[nodiscard]] bool isBoundary() const noexcept
    {
        return neighbor_ == nullptr;
    }

    [[nodiscard]] Polyhedron *ownerCell() noexcept
    {
        return owner_;
    }

    [[nodiscard]] const Polyhedron *
    ownerCell() const noexcept
    {
        return owner_;
    }

    [[nodiscard]] Polyhedron *
    neighborCell() noexcept
    {
        return neighbor_;
    }

    [[nodiscard]] const Polyhedron *
    neighborCell() const noexcept
    {
        return neighbor_;
    }

    [[nodiscard]] Face *neighborFace() noexcept
    {
        return neighborFace_;
    }

    [[nodiscard]] const Face *neighborFace() const noexcept
    {
        return neighborFace_;
    }

    void setOwnerCell(
        Polyhedron *cell) noexcept
    {
        owner_ = cell;
    }

    void setNeighborCell(
        Polyhedron *cell) noexcept
    {
        neighbor_ = cell;
    }

    void setNeighborFace(Face *face) noexcept
    {
        neighborFace_ = face;
    }

    [[nodiscard]] Point averagePoint() const;
    [[nodiscard]] Point centroid() const;
    [[nodiscard]] Point unitNormal() const;
    [[nodiscard]] double computeArea() const;

    /**
     * @brief 一次 face traversal 同时计算 area / centroid / normal。
     *
     * 供 CpGrid setup 使用；旧的三个独立查询接口继续保留。
     */
    [[nodiscard]] PolygonGeometry computeGeometry() const;

    void cacheArea(double value) noexcept
    {
        area_ = value;
    }

    void cacheCentroid(const Point &value) noexcept
    {
        centroid_ = value;
    }

    [[nodiscard]] const Point &cachedCentroid() const noexcept
    {
        return centroid_;
    }

    void cacheUnitNormal(const Point &value) noexcept
    {
        unitNormal_ = value;
    }

    [[nodiscard]] const Point &cachedUnitNormal() const noexcept
    {
        return unitNormal_;
    }

    [[nodiscard]] double area() const noexcept
    {
        return area_;
    }

    void cacheTransmissibility(
        double value) noexcept
    {
        transmissibility_ = value;
    }

    [[nodiscard]] double
    transmissibility() const noexcept
    {
        return transmissibility_;
    }

    void cacheGravityTerm(
        double value) noexcept
    {
        gravityTerm_ = value;
    }

    [[nodiscard]] double
    gravityTerm() const noexcept
    {
        return gravityTerm_;
    }

    friend std::ostream &operator<<(
        std::ostream &os,
        const Face &face);

  private:
    std::vector<Node *> nodes_;
    PetscInt inputPosition_{0};
    PetscInt inputFaceId_{0};

    Polyhedron *owner_{nullptr};
    Polyhedron *neighbor_{nullptr};
    Face *neighborFace_{nullptr};

    Point centroid_{};
    Point unitNormal_{};
    double area_{0.0};
    double transmissibility_{0.0};
    double gravityTerm_{0.0};
};

} // namespace MPMC
