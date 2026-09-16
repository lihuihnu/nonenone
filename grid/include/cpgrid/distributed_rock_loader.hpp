/**
 * @file distributed_rock_loader.hpp
 * @brief CpGrid root-only 岩石属性读取与分布式 PETSc Vec 装配。
 */
#pragma once

#include <cpgrid/cpgrid.hpp>
#include <cpgrid/csv_reader.hpp>
#include <cpgrid/grdecl.hpp>
#include <cpgrid/grdecl_rock_payload.hpp>
#include <cpgrid/rock_property_validation.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace MPMC
{

namespace detail
{

inline void broadcastRockRootError(
    MPI_Comm comm,
    int root,
    int rank,
    std::string &error)
{
    int success = rank == root && error.empty() ? 1 : 0;
    PetscCallMPIAbort(comm, MPI_Bcast(&success, 1, MPI_INT, root, comm));

    int length = rank == root ? static_cast<int>(error.size()) : 0;
    PetscCallMPIAbort(comm, MPI_Bcast(&length, 1, MPI_INT, root, comm));
    if (length < 0)
        throw std::runtime_error("Invalid root rock-input error length.");
    if (rank != root)
        error.resize(static_cast<std::size_t>(length));
    if (length > 0)
        PetscCallMPIAbort(comm, MPI_Bcast(error.data(), length, MPI_CHAR, root, comm));

    if (success == 0)
        throw std::runtime_error("Root-only rock ingest failed: " + error);
}

template <class Grid>
void applyRockPropertiesFromRootData(
    Grid &grid,
    const std::vector<double> &porosity,
    const std::array<std::vector<double>, 3> &permeability,
    const RockReplacementMeans &means,
    const CpGridRockLoadOptions &options,
    int root)
{
    if (grid.isSetup())
        throw std::logic_error("Rock properties cannot be replaced after CpGrid::setup().");

    const int rank = grid.mesh().rank();
    const DofMap &permMap = grid.dofMap(3);
    const DofMap &poroMap = grid.dofMap(1);
    Vec permVec = grid.permeabilityVector();
    Vec poroVec = grid.porosityVector();

    if (rank == root)
    {
        constexpr std::size_t kCellChunk = 4096;
        std::vector<PetscInt> poroIndices;
        std::vector<PetscScalar> poroValues;
        std::vector<PetscInt> permIndices;
        std::vector<PetscScalar> permValues;
        poroIndices.reserve(kCellChunk);
        poroValues.reserve(kCellChunk);
        permIndices.reserve(3 * kCellChunk);
        permValues.reserve(3 * kCellChunk);

        const auto flush = [&]
        {
            if (poroIndices.empty())
                return;
            PetscCallAbort(
                grid.mesh().communicator(),
                VecSetValues(
                    poroVec,
                    static_cast<PetscInt>(poroIndices.size()),
                    poroIndices.data(),
                    poroValues.data(),
                    INSERT_VALUES));
            PetscCallAbort(
                grid.mesh().communicator(),
                VecSetValues(
                    permVec,
                    static_cast<PetscInt>(permIndices.size()),
                    permIndices.data(),
                    permValues.data(),
                    INSERT_VALUES));
            poroIndices.clear();
            poroValues.clear();
            permIndices.clear();
            permValues.clear();
        };

        for (std::size_t input = 0; input < porosity.size(); ++input)
        {
            const PetscInt currentId = grid.mesh().currentIdFromInputIndex(input);
            double phi = porosity[input];
            if (phi == 0.0 && options.replaceZeroPorosityWithMean)
                phi = means.porosity;
            poroIndices.push_back(poroMap.globalIndex(currentId, 0));
            poroValues.push_back(static_cast<PetscScalar>(phi));

            for (PetscInt axis = 0; axis < 3; ++axis)
            {
                double k = permeability[static_cast<std::size_t>(axis)][input];
                if (k == 0.0 && options.replaceZeroPermeabilityWithMean)
                    k = means.permeability[static_cast<std::size_t>(axis)];
                permIndices.push_back(permMap.globalIndex(currentId, axis));
                permValues.push_back(static_cast<PetscScalar>(k));
            }
            if (poroIndices.size() == kCellChunk)
                flush();
        }
        flush();
    }

    PetscCallAbort(grid.mesh().communicator(), VecAssemblyBegin(permVec));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyEnd(permVec));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyBegin(poroVec));
    PetscCallAbort(grid.mesh().communicator(), VecAssemblyEnd(poroVec));
}

} // namespace detail

/**
 * @brief 仅 root 持有轻量岩石数组，并由 root 装配分布式 PETSc Vec。
 */
template <class Grid>
void loadRockPropertiesFromRootData(
    Grid &grid,
    const CpGridRootRockData *rockOnRoot,
    const CpGridRockLoadOptions &options = {},
    int root = 0)
{
    const int rank = grid.mesh().rank();
    const int size = grid.mesh().processCount();
    if (root < 0 || root >= size)
        throw std::out_of_range("Rock ingest root rank is outside communicator.");

    detail::RockReplacementMeans means;
    std::string error;

    if (rank == root)
    {
        try
        {
            if (rockOnRoot == nullptr)
                throw std::invalid_argument("Root rank requires rock property data.");
            if (rockOnRoot->cellCount != grid.mesh().cellCount())
                throw std::runtime_error(
                    "Rock property count must match active Mesh cell count.");
            means = detail::validateRockProperties(
                rockOnRoot->porosity,
                rockOnRoot->permeability,
                grid.mesh().cellCount(),
                options);
        }
        catch (const std::exception &caught)
        {
            error = caught.what();
        }
    }
    detail::broadcastRockRootError(grid.mesh().communicator(), root, rank, error);

    static const std::vector<double> emptyPorosity;
    static const std::array<std::vector<double>, 3> emptyPermeability{};
    detail::applyRockPropertiesFromRootData(
        grid,
        rank == root ? rockOnRoot->porosity : emptyPorosity,
        rank == root ? rockOnRoot->permeability : emptyPermeability,
        means,
        options,
        root);
}

/**
 * @brief 仅 root 读取 MRST rock CSV，并由 root 装配分布式 PETSc Vec。
 */
template <class Grid>
void loadRockPropertiesFromCsvRoot(
    Grid &grid,
    const std::string &dataDirectory,
    const CpGridRockLoadOptions &options = {},
    int root = 0)
{
    const int rank = grid.mesh().rank();
    const int size = grid.mesh().processCount();
    if (root < 0 || root >= size)
        throw std::out_of_range("Rock ingest root rank is outside communicator.");

    std::vector<double> porosity;
    std::array<std::vector<double>, 3> permeability;
    detail::RockReplacementMeans means;
    std::string error;

    if (rank == root)
    {
        try
        {
            std::string base = dataDirectory;
            while (!base.empty() && base.back() == '/')
                base.pop_back();
            const auto permTable = readCsv<double>(base + "/rock/perm/data.csv");
            const auto poroTable = readCsv<double>(base + "/rock/poro/data.csv");
            if (permTable.size() != grid.mesh().cellCount() ||
                poroTable.size() != grid.mesh().cellCount())
                throw std::runtime_error("Rock CSV row count must match active Mesh cell count.");

            porosity.resize(poroTable.size());
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
            means = detail::validateRockProperties(
                porosity, permeability, grid.mesh().cellCount(), options);
        }
        catch (const std::exception &caught)
        {
            error = caught.what();
        }
    }
    detail::broadcastRockRootError(grid.mesh().communicator(), root, rank, error);
    detail::applyRockPropertiesFromRootData(
        grid, porosity, permeability, means, options, root);
}

/**
 * @brief 仅 root 持有 GRDECL rock 数组，并由 root 装配分布式 PETSc Vec。
 */
template <class Grid>
void loadRockPropertiesFromGrdeclRoot(
    Grid &grid,
    const GrdeclGridData *grdeclOnRoot,
    const CpGridRockLoadOptions &options = {},
    int root = 0)
{
    const int rank = grid.mesh().rank();
    const int size = grid.mesh().processCount();
    if (root < 0 || root >= size)
        throw std::out_of_range("Rock ingest root rank is outside communicator.");

    std::vector<double> porosity;
    std::array<std::vector<double>, 3> permeability;
    detail::RockReplacementMeans means;
    std::string error;

    if (rank == root)
    {
        try
        {
            if (grdeclOnRoot == nullptr)
                throw std::invalid_argument("Root rank requires GRDECL rock data.");
            if (!grdeclOnRoot->hasPorosity())
                throw std::runtime_error("GRDECL rock loading requires PORO.");
            if (!grdeclOnRoot->hasPermeability())
                throw std::runtime_error(
                    "GRDECL rock loading requires PERMX (PERMY/PERMZ may fall back to PERMX).");

            porosity = grdeclOnRoot->porosity;
            permeability = {{
                grdeclOnRoot->permeabilityX,
                grdeclOnRoot->permeabilityY,
                grdeclOnRoot->permeabilityZ}};
            means = detail::validateRockProperties(
                porosity, permeability, grid.mesh().cellCount(), options);
        }
        catch (const std::exception &caught)
        {
            error = caught.what();
        }
    }
    detail::broadcastRockRootError(grid.mesh().communicator(), root, rank, error);
    detail::applyRockPropertiesFromRootData(
        grid, porosity, permeability, means, options, root);
}

} // namespace MPMC
