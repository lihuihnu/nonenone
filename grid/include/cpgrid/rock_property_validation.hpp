/**
 * @file rock_property_validation.hpp
 * @brief CpGrid 岩石属性计数、范围和零值替换均值的共享校验规则。
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace MPMC
{

/** Optional data-cleaning policy shared by CSV and GRDECL rock input. */
struct CpGridRockLoadOptions final
{
    bool replaceZeroPermeabilityWithMean{false};
    bool replaceZeroPorosityWithMean{false};
};

namespace detail
{

struct RockReplacementMeans final
{
    double porosity{0.0};
    std::array<double, 3> permeability{};
};

inline RockReplacementMeans validateRockProperties(
    const std::vector<double> &porosity,
    const std::array<std::vector<double>, 3> &permeability,
    std::size_t expectedCount,
    const CpGridRockLoadOptions &options)
{
    if (porosity.size() != expectedCount || permeability[0].size() != expectedCount ||
        permeability[1].size() != expectedCount || permeability[2].size() != expectedCount)
        throw std::runtime_error("Rock property count must match active Mesh cell count.");

    RockReplacementMeans means;
    std::size_t porosityCount = 0;
    std::array<std::size_t, 3> permeabilityCount{};

    for (std::size_t i = 0; i < expectedCount; ++i)
    {
        const double phi = porosity[i];
        if (!std::isfinite(phi) || phi < 0.0 || phi > 1.0)
            throw std::runtime_error("Porosity must be finite and inside [0,1].");
        if (phi > 0.0)
        {
            means.porosity += phi;
            ++porosityCount;
        }

        for (std::size_t axis = 0; axis < 3; ++axis)
        {
            const double k = permeability[axis][i];
            if (!std::isfinite(k) || k < 0.0)
                throw std::runtime_error("Permeability must be finite and non-negative.");
            if (k > 0.0)
            {
                means.permeability[axis] += k;
                ++permeabilityCount[axis];
            }
        }
    }

    if (options.replaceZeroPorosityWithMean)
    {
        if (porosityCount == 0)
            throw std::runtime_error("Cannot replace zero porosity: no positive sample exists.");
        means.porosity /= static_cast<double>(porosityCount);
    }
    if (options.replaceZeroPermeabilityWithMean)
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
            if (permeabilityCount[axis] == 0)
                throw std::runtime_error(
                    "Cannot replace zero permeability: an axis has no positive sample.");
            means.permeability[axis] /= static_cast<double>(permeabilityCount[axis]);
        }
    return means;
}

} // namespace detail
} // namespace MPMC
