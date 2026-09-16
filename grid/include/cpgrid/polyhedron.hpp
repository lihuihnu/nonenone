/**
 * @file polyhedron.hpp
 * @brief CpGrid 单元多面体及输入/当前编号信息。
 */
#pragma once

#include <cpgrid/face.hpp>
#include <cpgrid/point.hpp>

#include <petscsys.h>

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

namespace MPMC
{

/** @brief 单元体积与体积加权质心的一次遍历结果。 */
struct PolyhedronGeometry final
{
    double volume{0.0};
    Point centroid{};
};

/**
 * @brief Mesh 中的一个多面体单元。
 *
 * Polyhedron 拥有 Face 集合。由于 Face 内部保存指向 Mesh 拓扑对象的裸观察
 * 指针，Polyhedron 禁止复制。移动只允许发生在拓扑绑定之前的 Mesh 构建阶段。
 */
class Polyhedron final
{
  public:
    Polyhedron() = default;

    Polyhedron(
        std::vector<Face> faces,
        PetscInt storageIndex) noexcept
        : currentId_(storageIndex),
          storageIndex_(storageIndex),
          faces_(std::move(faces))
    {
        rebindFaceOwners();
    }

    Polyhedron(const Polyhedron &) = delete;
    Polyhedron &operator=(const Polyhedron &) = delete;

    Polyhedron(Polyhedron &&other) noexcept;
    Polyhedron &operator=(
        Polyhedron &&other) noexcept;

    ~Polyhedron() = default;

    [[nodiscard]] PetscInt id() const noexcept
    {
        return currentId_;
    }

    void setId(PetscInt value) noexcept
    {
        currentId_ = value;
    }

    /**
     * @brief 返回单元在原始网格数据文件中的行序号。
     *
     * MRST canonical input creates cells in the exact row order of
     * `G/cells/indexMap/data.csv`; GRDECL canonical input preserves active-cell
     * storage order. Therefore this id is invariant under METIS
     * partitioning and `resetCurrentIds()`.  It is the canonical external id
     * used by case files, well perforation tables and exported result vectors.
     */
    [[nodiscard]] PetscInt inputIndex() const noexcept
    {
        return storageIndex_;
    }

    /**
     * @brief 返回稳定内部存储位置；数值上与 inputIndex() 相同。
     *
     * Kept for mesh-internal topology code.  User-facing code should prefer
     * `inputIndex()` so the meaning of the number is explicit.
     */
    [[nodiscard]] PetscInt storageIndex() const noexcept
    {
        return storageIndex_;
    }

    [[nodiscard]] int processorId() const noexcept
    {
        return processorId_;
    }

    void setProcessorId(int value) noexcept
    {
        processorId_ = value;
    }

    [[nodiscard]] bool isOwnedBy(
        int rank) const noexcept
    {
        return processorId_ == rank;
    }

    [[nodiscard]] std::vector<Face> &
    faces() noexcept
    {
        return faces_;
    }

    [[nodiscard]] const std::vector<Face> &
    faces() const noexcept
    {
        return faces_;
    }

    [[nodiscard]] std::size_t
    faceCount() const noexcept
    {
        return faces_.size();
    }

    void rebindFaceOwners() noexcept;

    [[nodiscard]] Point vertexAverage() const;
    [[nodiscard]] Point centroid() const;
    [[nodiscard]] double computeVolume() const;

    /**
     * @brief 一次 cell traversal 同时计算 volume 与 centroid。
     *
     * 供 CpGrid setup 使用；旧的独立 computeVolume()/centroid() 接口保留。
     */
    [[nodiscard]] PolyhedronGeometry computeGeometry() const;

    void cacheVolume(double value) noexcept
    {
        volume_ = value;
    }

    void cacheCentroid(const Point &value) noexcept
    {
        centroid_ = value;
    }

    [[nodiscard]] const Point &cachedCentroid() const noexcept
    {
        return centroid_;
    }

    [[nodiscard]] double volume() const noexcept
    {
        return volume_;
    }

    friend std::ostream &operator<<(
        std::ostream &os,
        const Polyhedron &cell);

  private:
    void resetMovedFrom_() noexcept;

    PetscInt currentId_{0};
    PetscInt storageIndex_{0};
    std::vector<Face> faces_;
    int processorId_{0};
    Point centroid_{};
    double volume_{0.0};
};

} // namespace MPMC
