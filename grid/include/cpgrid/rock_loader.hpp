/**
 * @file rock_loader.hpp
 * @brief 从标准化 CSV 读取孔隙度和渗透率等岩石属性。
 */
#pragma once

#include <cpgrid/cpgrid.hpp>
#include <cpgrid/csv_reader.hpp>
#include <cpgrid/grdecl.hpp>
#include <cpgrid/rock_property_validation.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace MPMC
{

namespace detail
{
template <class Grid>
void applyRockProperties(
    Grid &grid,
    const std::vector<double> &porosity,
    const std::array<std::vector<double>, 3> &permeability,
    const CpGridRockLoadOptions &options)
{
    if (grid.isSetup())
        throw std::logic_error("Rock properties cannot be replaced after CpGrid::setup().");

    const std::size_t n = grid.mesh().cellCount();
    const RockReplacementMeans means =
        validateRockProperties(porosity, permeability, n, options);

    const DofMap &permMap = grid.dofMap(3);
    const DofMap &poroMap = grid.dofMap(1);
    std::vector<PetscInt> permIndices;
    std::vector<PetscScalar> permValues;
    std::vector<PetscInt> poroIndices;
    std::vector<PetscScalar> poroValues;
    permIndices.reserve(static_cast<std::size_t>(grid.mesh().localCellCount()) * 3);
    permValues.reserve(permIndices.capacity());
    poroIndices.reserve(static_cast<std::size_t>(grid.mesh().localCellCount()));
    poroValues.reserve(poroIndices.capacity());

    for (PetscInt currentId : grid.mesh().ownedCellIds())
    {
        const Polyhedron &cell = grid.mesh().cellByCurrentId(currentId);
        const std::size_t i = static_cast<std::size_t>(cell.storageIndex());
        double phi = porosity[i];
        if (phi == 0.0 && options.replaceZeroPorosityWithMean)
            phi = means.porosity;
        poroIndices.push_back(poroMap.globalIndex(cell, 0));
        poroValues.push_back(static_cast<PetscScalar>(phi));

        for (PetscInt axis = 0; axis < 3; ++axis)
        {
            double k = permeability[static_cast<std::size_t>(axis)][i];
            if (k == 0.0 && options.replaceZeroPermeabilityWithMean)
                k = means.permeability[static_cast<std::size_t>(axis)];
            permIndices.push_back(permMap.globalIndex(cell, axis));
            permValues.push_back(static_cast<PetscScalar>(k));
        }
    }

    Vec permVec = grid.permeabilityVector();
    Vec poroVec = grid.porosityVector();
    PetscCallAbort(
        grid.mesh().communicator(),
        VecSetValues(
            permVec, static_cast<PetscInt>(permIndices.size()),
            permIndices.data(), permValues.data(), INSERT_VALUES));
    PetscCallAbort(
        grid.mesh().communicator(),
        VecSetValues(
            poroVec, static_cast<PetscInt>(poroIndices.size()),
            poroIndices.data(), poroValues.data(), INSERT_VALUES));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyBegin(permVec));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyEnd(permVec));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyBegin(poroVec));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyEnd(poroVec));
}
} // namespace detail

/** Load MRST-exported rock CSVs. Units are passed through unchanged. */
template <class Grid>
void loadRockPropertiesFromCsv(
    Grid &grid,
    const std::string &dataDirectory,
    const CpGridRockLoadOptions &options = {})
{
    std::string base = dataDirectory;
    while (!base.empty() && base.back() == '/')
        base.pop_back();

    const auto permTable = readCsv<double>(base + "/rock/perm/data.csv");
    const auto poroTable = readCsv<double>(base + "/rock/poro/data.csv");
    if (permTable.size() != grid.mesh().cellCount() ||
        poroTable.size() != grid.mesh().cellCount())
        throw std::runtime_error("Rock CSV row count must match active Mesh cell count.");

    std::vector<double> porosity(poroTable.size());
    std::array<std::vector<double>, 3> permeability;
    for (auto &values : permeability)
        values.resize(permTable.size());

    for (std::size_t i = 0; i < permTable.size(); ++i)
    {
        if (permTable[i].size() < 3 || poroTable[i].empty())
            throw std::runtime_error("Rock CSV contains an incomplete row.");
        porosity[i] = poroTable[i][0];
        for (std::size_t axis = 0; axis < 3; ++axis)
            permeability[axis][i] = permTable[i][axis];
    }

    detail::applyRockProperties(grid, porosity, permeability, options);
}

/** Load PORO/PERMX/PERMY/PERMZ parsed directly from a GRDECL file. */
template <class Grid>
void loadRockPropertiesFromGrdecl(
    Grid &grid,
    const GrdeclGridData &grdecl,
    const CpGridRockLoadOptions &options = {})
{
    if (!grdecl.hasPorosity())
        throw std::runtime_error("GRDECL rock loading requires PORO.");
    if (!grdecl.hasPermeability())
        throw std::runtime_error(
            "GRDECL rock loading requires PERMX (PERMY/PERMZ may fall back to PERMX).");

    const std::array<std::vector<double>, 3> permeability{{
        grdecl.permeabilityX,
        grdecl.permeabilityY,
        grdecl.permeabilityZ}};
    detail::applyRockProperties(grid, grdecl.porosity, permeability, options);
}

} // namespace MPMC
