/**
 * @file distributed_mesh_test.cpp
 * @brief 集成测试：验证 root-only CpGrid ingest 的 MPI owned+ghost 拓扑分发。
 */
#include <cpgrid/grdecl.hpp>
#include <cpgrid/mesh.hpp>

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

MPMC::GrdeclCell makeCell(
    int cartesianId,
    double x0,
    double x1,
    int left,
    int right)
{
    MPMC::GrdeclCell cell;
    cell.cartesianId = cartesianId;
    cell.neighborStorage = {{left, right, -1, -1, -1, -1}};
    for (int iz = 0; iz < 2; ++iz)
        for (int iy = 0; iy < 2; ++iy)
            for (int ix = 0; ix < 2; ++ix)
            {
                const int corner = ix + 2 * iy + 4 * iz;
                cell.corner[static_cast<std::size_t>(corner)] = {
                    ix ? x1 : x0,
                    static_cast<double>(iy),
                    static_cast<double>(iz)};
            }
    return cell;
}

MPMC::GrdeclGridData makeLinearGrid(int cellCount)
{
    MPMC::GrdeclGridData grid;
    grid.sourcePath = "synthetic-root-only.grdecl";
    grid.nx = cellCount;
    grid.ny = 1;
    grid.nz = 1;
    grid.nodeMergeTolerance = 1.0e-12;
    grid.activeCells.reserve(static_cast<std::size_t>(cellCount));
    for (int cell = 0; cell < cellCount; ++cell)
        grid.activeCells.push_back(
            makeCell(
                cell,
                static_cast<double>(cell),
                static_cast<double>(cell + 1),
                cell == 0 ? -1 : cell - 1,
                cell + 1 == cellCount ? -1 : cell + 1));
    return grid;
}

} // namespace

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    const int globalCells = std::max(8, 4 * size);
    std::unique_ptr<MPMC::Mesh> rootMesh;
    if (rank == 0)
    {
        auto input = makeLinearGrid(globalCells);
        rootMesh = std::make_unique<MPMC::Mesh>(input, PETSC_COMM_WORLD);
        rootMesh->prepareForRootDistribution(0);
    }

    auto mesh = MPMC::Mesh::distributePreparedRoot(
        std::move(rootMesh), PETSC_COMM_WORLD, 0);

    double localError = 0.0;
    if (mesh->cellCount() != static_cast<std::size_t>(globalCells))
        localError = 1.0;
    if (rank == 0 && mesh->isLocalSnapshotCompacted())
        localError = 1.0;
    if (rank != 0 && !mesh->isLocalSnapshotCompacted())
        localError = 1.0;

    // 生产 DQ 路径在一次性全局元数据写出后也会释放 root 的远端拓扑。
    mesh->compactToLocalSnapshot(false);
    if (size > 1 && !mesh->isLocalSnapshotCompacted())
        localError = 1.0;
    if (mesh->materializedCellCount() !=
        mesh->ownedCellIds().size() + mesh->ghostCellIds().size())
        localError = 1.0;

    std::vector<PetscInt> seen(static_cast<std::size_t>(globalCells), PetscInt(-1));
    for (int input = 0; input < globalCells; ++input)
    {
        const PetscInt current = mesh->currentIdFromInputIndex(static_cast<std::size_t>(input));
        if (current < 0 || current >= globalCells)
        {
            localError = 1.0;
            continue;
        }
        if (seen[static_cast<std::size_t>(current)] >= 0)
            localError = 1.0;
        seen[static_cast<std::size_t>(current)] = input;
        if (mesh->inputIndexFromCurrentId(current) != input)
            localError = 1.0;
        if (mesh->currentIdFromCartesianId(input) != current)
            localError = 1.0;
    }

    for (PetscInt current : mesh->ownedCellIds())
    {
        const auto &cell = mesh->cellByCurrentId(current);
        if (!cell.isOwnedBy(rank))
            localError = 1.0;
        if (mesh->ownerRankFromCurrentId(current) != rank)
            localError = 1.0;
        localError = std::max(localError, std::abs(cell.computeVolume() - 1.0));

        for (const auto &face : cell.faces())
        {
            const auto *neighbor = face.neighborCell();
            if (neighbor == nullptr)
                continue;
            if (mesh->ownerRankFromCurrentId(neighbor->id()) != neighbor->processorId())
                localError = 1.0;
            try
            {
                (void)mesh->cellByCurrentId(neighbor->id());
            }
            catch (...)
            {
                localError = 1.0;
            }
            if (face.neighborFace() == nullptr)
                localError = 1.0;
        }
    }

    double globalError = 0.0;
    MPI_Allreduce(
        &localError,
        &globalError,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        PETSC_COMM_WORLD);

    if (rank == 0)
    {
        PetscPrintf(
            PETSC_COMM_SELF,
            "root-only distributed Mesh: cells=%d ranks=%d status=%s\n",
            globalCells,
            size,
            globalError <= 1.0e-12 ? "PASS" : "FAIL");
    }

    PetscFinalize();
    return globalError <= 1.0e-12 ? 0 : 1;
}
