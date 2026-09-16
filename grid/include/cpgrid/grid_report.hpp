/**
 * @file grid_report.hpp
 * @brief 网格规模、分区和几何信息的统一报告接口。
 */
#pragma once

#include <common/console.hpp>
#include <cpgrid/cpgrid.hpp>
#include <grid/report_statistics.hpp>

#include <petscsys.h>

#include <array>
#include <sstream>
#include <stdexcept>
#include <string>

namespace MPMC
{

/**
 * Collective CpGrid report. All ranks enter; only rank 0 prints.
 * Rock properties and cell volumes are summarized over owned cells.
 */
template <class Grid>
void printCpGridSummary(const Grid &grid)
{
    const Mesh &mesh = grid.mesh();
    const auto rock = grid.rockLocalView();

    constexpr double squareMetreToMilliDarcy = 1.0 / 9.869232667160130e-16;
    detail::GridReportStatistics localStatistics;

    for (PetscInt currentId : mesh.ownedCellIds())
    {
        const auto &cell = mesh.cellByCurrentId(currentId);
        const double volume = grid.cellVolume(cell);
        const double porosity = rock.porosity(cell);
        const auto permeability = rock.permeability(cell);
        localStatistics.add(volume, porosity, permeability);
    }

    const detail::GridReportStatistics globalStatistics =
        detail::reduceGridReportStatistics(localStatistics, mesh.communicator());

    if (mesh.rank() != 0)
        return;
    if (globalStatistics.count <= 0)
        throw std::runtime_error("CpGrid report found no owned cells.");

    const auto dims = mesh.logicalDimensions();
    ConsoleSection section("GRID LOAD / ROCK PROPERTIES");
    section.row("Load status", "SUCCESS")
        .row("Input format", mesh.sourceFormat())
        .row("Data location", mesh.sourcePath())
        .row("MPI ranks", std::to_string(mesh.processCount()));
    if (dims[0] > 0 && dims[1] > 0 && dims[2] > 0)
    {
        std::ostringstream logical;
        logical << dims[0] << " x " << dims[1] << " x " << dims[2];
        section.row("Logical size", logical.str());
    }
    section.row("Active cells", std::to_string(mesh.cellCount()))
        .row("Nodes", std::to_string(mesh.nodeCount()))
        .row("Unique faces", std::to_string(mesh.uniqueFaceCount()))
        .row("Topology", mesh.isTopologyInitialized() ? "READY" : "NOT READY")
        .separator()
        .row("Bulk volume", consoleNumber(globalStatistics.volumeSum), "m3")
        .row("Cell volume", formatGridStatistics(
            globalStatistics.volumeStatistics()), "m3")
        .row("Porosity", formatGridStatistics(
            globalStatistics.porosityStatistics()));

    static constexpr std::array<const char *, 3> axisName{{"Kx", "Ky", "Kz"}};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        section.row(
            axisName[axis],
            formatGridStatistics({
                globalStatistics.permeabilityMin[axis] * squareMetreToMilliDarcy,
                globalStatistics.permeabilityStatistics(axis).average * squareMetreToMilliDarcy,
                globalStatistics.permeabilityMax[axis] * squareMetreToMilliDarcy}),
            "mD");
    }

    const std::string text = section.str();
    PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
}

} // namespace MPMC
