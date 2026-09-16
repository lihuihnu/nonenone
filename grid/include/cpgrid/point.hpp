/**
 * @file point.hpp
 * @brief 三维点和向量几何基础类型。
 */
#pragma once

#include <array>
#include <cstddef>
#include <iosfwd>

namespace MPMC
{

/**
 * @brief 三维笛卡尔几何向量。
 *
 * 即使 MESH_DIM=2，也保留 z 分量并令其通常为 0。
 * Point 只表示几何量，不承担网格拓扑或物理属性职责。
 */
class Point
{
  public:
    constexpr Point() noexcept = default;

    constexpr Point(
        double x,
        double y,
        double z = 0.0) noexcept
        : coords_{x, y, z}
    {
    }

    [[nodiscard]] constexpr double x() const noexcept
    {
        return coords_[0];
    }

    [[nodiscard]] constexpr double y() const noexcept
    {
        return coords_[1];
    }

    [[nodiscard]] constexpr double z() const noexcept
    {
        return coords_[2];
    }

    [[nodiscard]] constexpr double operator[](
        std::size_t index) const noexcept
    {
        return coords_[index];
    }

    [[nodiscard]] constexpr double &operator[](
        std::size_t index) noexcept
    {
        return coords_[index];
    }

    [[nodiscard]] Point operator+(
        const Point &rhs) const noexcept;

    [[nodiscard]] Point operator-(
        const Point &rhs) const noexcept;

    [[nodiscard]] Point operator*(
        double scalar) const noexcept;

    [[nodiscard]] Point operator/(
        double scalar) const;

    Point &operator+=(const Point &rhs) noexcept;
    Point &operator-=(const Point &rhs) noexcept;
    Point &operator*=(double scalar) noexcept;
    Point &operator/=(double scalar);

    [[nodiscard]] double dot(
        const Point &rhs) const noexcept;

    [[nodiscard]] double normSquared() const noexcept;
    [[nodiscard]] double norm() const noexcept;

    [[nodiscard]] bool isFinite() const noexcept;

  private:
    std::array<double, 3> coords_{0.0, 0.0, 0.0};
};

[[nodiscard]] Point operator*(
    double scalar,
    const Point &point) noexcept;

std::ostream &operator<<(
    std::ostream &os,
    const Point &point);

} // namespace MPMC
