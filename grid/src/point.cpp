/**
 * @file point.cpp
 * @brief 三维点和向量几何基础类型的实现。
 */
#include <cpgrid/point.hpp>

#include <cmath>
#include <ostream>
#include <stdexcept>

namespace MPMC
{

Point Point::operator+(
    const Point &rhs) const noexcept
{
    return Point(
        x() + rhs.x(),
        y() + rhs.y(),
        z() + rhs.z());
}

Point Point::operator-(
    const Point &rhs) const noexcept
{
    return Point(
        x() - rhs.x(),
        y() - rhs.y(),
        z() - rhs.z());
}

Point Point::operator*(
    double scalar) const noexcept
{
    return Point(
        x() * scalar,
        y() * scalar,
        z() * scalar);
}

Point Point::operator/(
    double scalar) const
{
    if (scalar == 0.0)
    {
        throw std::invalid_argument(
            "Point division by zero.");
    }

    return Point(
        x() / scalar,
        y() / scalar,
        z() / scalar);
}

Point &Point::operator+=(
    const Point &rhs) noexcept
{
    coords_[0] += rhs.x();
    coords_[1] += rhs.y();
    coords_[2] += rhs.z();
    return *this;
}

Point &Point::operator-=(
    const Point &rhs) noexcept
{
    coords_[0] -= rhs.x();
    coords_[1] -= rhs.y();
    coords_[2] -= rhs.z();
    return *this;
}

Point &Point::operator*=(
    double scalar) noexcept
{
    coords_[0] *= scalar;
    coords_[1] *= scalar;
    coords_[2] *= scalar;
    return *this;
}

Point &Point::operator/=(
    double scalar)
{
    if (scalar == 0.0)
    {
        throw std::invalid_argument(
            "Point division by zero.");
    }

    coords_[0] /= scalar;
    coords_[1] /= scalar;
    coords_[2] /= scalar;
    return *this;
}

double Point::dot(
    const Point &rhs) const noexcept
{
    return x() * rhs.x() +
           y() * rhs.y() +
           z() * rhs.z();
}

double Point::normSquared() const noexcept
{
    return dot(*this);
}

double Point::norm() const noexcept
{
    return std::sqrt(normSquared());
}

bool Point::isFinite() const noexcept
{
    return std::isfinite(x()) && std::isfinite(y()) && std::isfinite(z());
}

Point operator*(
    double scalar,
    const Point &point) noexcept
{
    return point * scalar;
}

std::ostream &operator<<(
    std::ostream &os,
    const Point &point)
{
    return os
           << '('
           << point.x() << ", "
           << point.y() << ", "
           << point.z()
           << ')';
}

} // namespace MPMC
