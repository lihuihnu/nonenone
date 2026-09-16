/**
 * @file report_statistics.hpp
 * @brief 网格报告共享的统计累积、MPI 归约与格式化辅助。
 */
#pragma once

#include <common/console.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>

namespace MPMC
{

struct GridScalarStatistics final
{
    double minimum{0.0};
    double average{0.0};
    double maximum{0.0};
};

inline std::string formatGridStatistics(
    const GridScalarStatistics &stats,
    std::string_view unit = {})
{
    std::ostringstream out;
    out << "min=" << consoleNumber(stats.minimum)
        << "  avg=" << consoleNumber(stats.average)
        << "  max=" << consoleNumber(stats.maximum);
    if (!unit.empty())
        out << ' ' << unit;
    return out.str();
}

namespace detail
{

struct GridReportStatistics final
{
    long long count{0};
    double volumeMin{std::numeric_limits<double>::infinity()};
    double volumeMax{-std::numeric_limits<double>::infinity()};
    double volumeSum{0.0};
    double porosityMin{std::numeric_limits<double>::infinity()};
    double porosityMax{-std::numeric_limits<double>::infinity()};
    double porositySum{0.0};
    std::array<double, 3> permeabilityMin{{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()}};
    std::array<double, 3> permeabilityMax{{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()}};
    std::array<double, 3> permeabilitySum{};

    void add(double volume, double porosity, const std::array<double, 3> &permeability) noexcept
    {
        volumeMin = std::min(volumeMin, volume);
        volumeMax = std::max(volumeMax, volume);
        volumeSum += volume;
        porosityMin = std::min(porosityMin, porosity);
        porosityMax = std::max(porosityMax, porosity);
        porositySum += porosity;
        for (std::size_t axis = 0; axis < permeability.size(); ++axis)
        {
            permeabilityMin[axis] = std::min(permeabilityMin[axis], permeability[axis]);
            permeabilityMax[axis] = std::max(permeabilityMax[axis], permeability[axis]);
            permeabilitySum[axis] += permeability[axis];
        }
        ++count;
    }

    [[nodiscard]] GridScalarStatistics volumeStatistics() const noexcept
    {
        return {volumeMin, volumeSum / static_cast<double>(count), volumeMax};
    }

    [[nodiscard]] GridScalarStatistics porosityStatistics() const noexcept
    {
        return {porosityMin, porositySum / static_cast<double>(count), porosityMax};
    }

    [[nodiscard]] GridScalarStatistics permeabilityStatistics(std::size_t axis) const noexcept
    {
        return {permeabilityMin[axis],
                permeabilitySum[axis] / static_cast<double>(count),
                permeabilityMax[axis]};
    }
};

inline GridReportStatistics reduceGridReportStatistics(
    const GridReportStatistics &local,
    MPI_Comm communicator)
{
    GridReportStatistics global;
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.count, &global.count, 1, MPI_LONG_LONG_INT, MPI_SUM, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.volumeMin, &global.volumeMin, 1, MPI_DOUBLE, MPI_MIN, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.volumeMax, &global.volumeMax, 1, MPI_DOUBLE, MPI_MAX, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.volumeSum, &global.volumeSum, 1, MPI_DOUBLE, MPI_SUM, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.porosityMin, &global.porosityMin, 1, MPI_DOUBLE, MPI_MIN, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.porosityMax, &global.porosityMax, 1, MPI_DOUBLE, MPI_MAX, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        &local.porositySum, &global.porositySum, 1, MPI_DOUBLE, MPI_SUM, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        local.permeabilityMin.data(), global.permeabilityMin.data(),
        3, MPI_DOUBLE, MPI_MIN, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        local.permeabilityMax.data(), global.permeabilityMax.data(),
        3, MPI_DOUBLE, MPI_MAX, communicator));
    PetscCallMPIAbort(communicator, MPI_Allreduce(
        local.permeabilitySum.data(), global.permeabilitySum.data(),
        3, MPI_DOUBLE, MPI_SUM, communicator));
    return global;
}

} // namespace detail
} // namespace MPMC
