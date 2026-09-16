/**
 * @file sparsity.cpp
 * @brief 由网格邻接关系构造 PETSc Jacobian 稀疏结构的实现。
 */
#include <cpgrid/sparsity.hpp>

#include <cpgrid/face.hpp>
#include <cpgrid/polyhedron.hpp>
#include <cpgrid/mesh.hpp>

#include <petscmat.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace MPMC
{

Sparsity::Sparsity(
    const DofMap &dofMap)
    : dofMap_(dofMap)
{
    build_();
}

void Sparsity::build_()
{
    const Mesh &mesh = dofMap_.mesh();
    const int rank = mesh.rank();

    diagonalBlockNnz_.clear();
    offDiagonalBlockNnz_.clear();

    diagonalBlockNnz_.reserve(
        static_cast<std::size_t>(
            mesh.localCellCount()));
    offDiagonalBlockNnz_.reserve(
        static_cast<std::size_t>(
            mesh.localCellCount()));

    // Reservoir cells have only a handful of faces. Reusing two compact
    // vectors avoids constructing two unordered-set hash tables per cell
    // while still protecting preallocation from duplicate face entries.
    std::vector<PetscInt> diagonalBlocks;
    std::vector<PetscInt> offDiagonalBlocks;
    diagonalBlocks.reserve(8);
    offDiagonalBlocks.reserve(8);

    const auto insertUnique = [](std::vector<PetscInt> &blocks, PetscInt id) {
        if (std::find(blocks.begin(), blocks.end(), id) == blocks.end())
            blocks.push_back(id);
    };

    for (PetscInt cellId : mesh.ownedCellIds())
    {
        const Polyhedron &cell = mesh.cellByCurrentId(cellId);
        diagonalBlocks.clear();
        offDiagonalBlocks.clear();

        diagonalBlocks.push_back(cell.id());

        for (const Face &face :
             cell.faces())
        {
            const Polyhedron *neighbor =
                face.neighborCell();

            if (neighbor == nullptr)
            {
                continue;
            }

            if (neighbor->isOwnedBy(rank))
            {
                insertUnique(diagonalBlocks, neighbor->id());
            }
            else
            {
                insertUnique(offDiagonalBlocks, neighbor->id());
            }
        }

        diagonalBlockNnz_.push_back(
            static_cast<PetscInt>(
                diagonalBlocks.size()));

        offDiagonalBlockNnz_.push_back(
            static_cast<PetscInt>(
                offDiagonalBlocks.size()));
    }

    if (diagonalBlockNnz_.size() !=
        static_cast<std::size_t>(
            mesh.localCellCount()))
    {
        throw std::runtime_error(
            "Sparsity did not build one block row per local cell.");
    }
}

PetscErrorCode Sparsity::createMatrix(
    Mat *matrix) const
{
    PetscFunctionBegin;

    const Mesh &mesh = dofMap_.mesh();

    PetscCheck(
        matrix != nullptr,
        mesh.communicator(),
        PETSC_ERR_ARG_NULL,
        "Output Mat pointer must not be null.");

    *matrix = nullptr;

    PetscCall(
        MatCreate(
            mesh.communicator(),
            matrix));

    PetscCall(
        MatSetSizes(
            *matrix,
            dofMap_.ownedDofCount(),
            dofMap_.ownedDofCount(),
            dofMap_.globalDofCount(),
            dofMap_.globalDofCount()));

    PetscCall(
        MatSetBlockSize(
            *matrix,
            dofMap_.dofPerCell()));

    /*
     * BAIJ 与当前“一个 cell 一个 dense block”的离散最匹配。
     * MatSetFromOptions() 允许用户用 -mat_type aij 等覆盖默认类型。
     */
    PetscCall(
        MatSetType(
            *matrix,
            MATBAIJ));

    PetscCall(
        MatSetFromOptions(
            *matrix));

    PetscCall(
        MatXAIJSetPreallocation(
            *matrix,
            dofMap_.dofPerCell(),
            diagonalBlockNnz_.empty()
                ? nullptr
                : diagonalBlockNnz_.data(),
            offDiagonalBlockNnz_.empty()
                ? nullptr
                : offDiagonalBlockNnz_.data(),
            nullptr,
            nullptr));

    ISLocalToGlobalMapping mapping =
        dofMap_.createLocalToGlobalMapping();

    PetscCall(
        MatSetLocalToGlobalMapping(
            *matrix,
            mapping,
            mapping));

    PetscCall(
        ISLocalToGlobalMappingDestroy(
            &mapping));

    PetscCall(
        MatSetOption(
            *matrix,
            MAT_NEW_NONZERO_ALLOCATION_ERR,
            PETSC_TRUE));

    PetscCall(
        MatSetUp(*matrix));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace MPMC
