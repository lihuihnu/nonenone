/**
 * @file grdecl_geometry.cpp
 * @brief GRDECL 长度/渗透率单位、MAPAXES 与 SI 标准化实现。
 */
#include "grdecl_detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace MPMC::grdecl_detail
{
namespace
{
constexpr double milliDarcyToSquareMetre = 9.869232667160130e-16;
constexpr double footToMetre = 0.3048;
constexpr double centimetreToMetre = 0.01;
} // namespace

double geometryLengthFactor(GrdeclUnitSystem units) noexcept
{
    if (units == GrdeclUnitSystem::Field)
        return footToMetre;
    if (units == GrdeclUnitSystem::Lab)
        return centimetreToMetre;
    return 1.0;
}

void convertGeometryToSi(
    std::vector<double> &coord,
    std::vector<double> &zcorn,
    GrdeclUnitSystem units)
{
    const double factor = geometryLengthFactor(units);
    for (double &value : coord)
        value *= factor;
    for (double &value : zcorn)
        value *= factor;
}

void applyMapAxes(
    std::vector<double> &coord,
    const std::vector<double> &mapAxes,
    double originLengthFactor)
{
    if (mapAxes.size() != 6)
        throw std::runtime_error("MAPAXES must contain exactly six values: X1 Y1 X2 Y2 X3 Y3.");
    if (!std::all_of(mapAxes.begin(), mapAxes.end(), [](double value) {
            return std::isfinite(value);
        }))
        throw std::runtime_error("MAPAXES values must be finite.");

    // 数值：X2/Y2 为原点，X2->X3 与 X2->X1 分别定义局部 X/Y 方向。
    double unitXx = mapAxes[4] - mapAxes[2];
    double unitXy = mapAxes[5] - mapAxes[3];
    double unitYx = mapAxes[0] - mapAxes[2];
    double unitYy = mapAxes[1] - mapAxes[3];
    const double normX = std::hypot(unitXx, unitXy);
    const double normY = std::hypot(unitYx, unitYy);
    if (!(normX > 0.0) || !(normY > 0.0) ||
        !std::isfinite(normX) || !std::isfinite(normY))
        throw std::runtime_error("MAPAXES must define two non-zero axes.");
    unitXx /= normX;
    unitXy /= normX;
    unitYx /= normY;
    unitYy /= normY;
    const double determinant = unitXx * unitYy - unitXy * unitYx;
    if (!std::isfinite(determinant) || std::abs(determinant) <= 1.0e-12)
        throw std::runtime_error("MAPAXES axes must not be collinear.");

    const double originX = originLengthFactor * mapAxes[2];
    const double originY = originLengthFactor * mapAxes[3];
    for (std::size_t index = 0; index < coord.size(); index += 3)
    {
        const double localX = coord[index];
        const double localY = coord[index + 1];
        coord[index] = originX + localX * unitXx + localY * unitYx;
        coord[index + 1] = originY + localX * unitXy + localY * unitYy;
    }
}

void convertPermeabilityToSi(std::vector<double> &values) noexcept
{
    for (double &value : values)
        value *= milliDarcyToSquareMetre;
}

} // namespace MPMC::grdecl_detail
