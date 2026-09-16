/**
 * @file grid_report.hpp
 * @brief 网格规模、分区和几何信息的统一报告接口。
 */
#pragma once

#include <common/console.hpp>
#include <grid/report_statistics.hpp>
#include <structuredgrid/structuredgrid.hpp>

#include <petscsys.h>

#include <array>
#include <sstream>
#include <stdexcept>
#include <string>

namespace MPMC
{

/** Collective report for generated StructuredGrid meshes. */
template <class Grid>
void printStructuredGridSummary(const Grid &grid)
{
    if (!grid.isSetup())
        throw std::logic_error("StructuredGrid report requires grid.setup().");
    if (!grid.hasRockProperties())
        throw std::logic_error("StructuredGrid report requires initialized rock properties.");

    const MPI_Comm comm = grid.communicator();

    constexpr double squareMetreToMilliDarcy = 1.0 / 9.869232667160130e-16;
    detail::GridReportStatistics localStatistics;
    {
        auto rock = grid.rockReadView();
        for (int global : grid.cellIndices())
        {
            const auto ijk = grid.globalToIJK(global);
            const int i = ijk[0];
            const int j = ijk[1];
            const int k = ijk[2];
            const auto cellSize = grid.cellSize(ijk);
            const auto permeability = rock.permeability(i, j, k);
            localStatistics.add(
                cellSize[0] * cellSize[1] * cellSize[2],
                rock.porosity(i, j, k),
                permeability);
        }
    }

    const detail::GridReportStatistics globalStatistics =
        detail::reduceGridReportStatistics(localStatistics, comm);

    PetscMPIInt rank = 0;
    PetscMPIInt size = 1;
    PetscCallMPIAbort(comm, MPI_Comm_rank(comm, &rank));
    PetscCallMPIAbort(comm, MPI_Comm_size(comm, &size));
    if (rank != 0)
        return;
    if (globalStatistics.count <= 0)
        throw std::runtime_error("StructuredGrid report found no active cells.");

    const auto dims = grid.dimensions();
    std::ostringstream logical;
    logical << dims[0] << " x " << dims[1] << " x " << dims[2];
    const auto physicalExtent = grid.physicalExtent();
    std::ostringstream extent;
    extent << consoleNumber(physicalExtent[0]) << " x "
           << consoleNumber(physicalExtent[1]) << " x "
           << consoleNumber(physicalExtent[2]) << " m";

    ConsoleSection section("GRID LOAD / ROCK PROPERTIES");
    section.row("Load status", "SUCCESS")
        .row("Input format", "StructuredGrid (generated)")
        .row("Data location", "case_config.hpp / in-memory")
        .row("MPI ranks", std::to_string(size))
        .row("Logical size", logical.str())
        .row("Physical extent", extent.str())
        .row("Active cells", std::to_string(globalStatistics.count))
        .row("Topology", "READY (DMDA Cartesian)")
        .separator()
        .row("Bulk volume", consoleNumber(globalStatistics.volumeSum), "m3")
        .row("Cell volume", formatGridStatistics(
            globalStatistics.volumeStatistics()), "m3")
        .row("Porosity", formatGridStatistics(
            globalStatistics.porosityStatistics()));

    static constexpr std::array<const char *, 3> axisName{{"Kx", "Ky", "Kz"}};
    for (std::size_t axis = 0; axis < 3; ++axis)
        section.row(
            axisName[axis],
            formatGridStatistics({
                globalStatistics.permeabilityMin[axis] * squareMetreToMilliDarcy,
                globalStatistics.permeabilityStatistics(axis).average * squareMetreToMilliDarcy,
                globalStatistics.permeabilityMax[axis] * squareMetreToMilliDarcy}),
            "mD");

    const std::string text = section.str();
    PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
}

} // namespace MPMC
