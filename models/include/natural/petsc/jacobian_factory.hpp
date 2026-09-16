/**
 * @file jacobian_factory.hpp
 * @brief 按静态网格/井拓扑创建、预分配并物化 PETSc Jacobian 结构。
 */
#pragma once

#include <natural/petsc/well_runtime.hpp>

#include <petscmat.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace MPMC
{

/**
 * @brief Natural Jacobian 的静态 block-CSR 拓扑。
 *
 * rows 按当前 rank 的 PETSc owned block ordinal 排列；columnBlocks 保存全局
 * block id。该结构只依赖 mesh partition、DOF layout 与井 perforation topology，
 * 与 Newton state、井开关和 RATE/BHP 当前控制状态无关。
 */
struct NaturalJacobianPattern final
{
    PetscInt blockSize{0};
    std::vector<PetscInt> rowBlocks;
    std::vector<PetscInt> rowOffsets;
    std::vector<PetscInt> columnBlocks;
    std::vector<PetscInt> diagonalNnz;
    std::vector<PetscInt> offDiagonalNnz;
};

/** @brief setup 阶段一次构造 Natural Jacobian 的精确 block topology。 */
template <class Indices, class Backend>
[[nodiscard]] NaturalJacobianPattern buildNaturalJacobianPattern(
    const Backend &backend,
    const std::vector<NaturalWell<Indices>> &wells)
{
    static_assert(
        Indices::numPrimaryVariables == Indices::numEquations,
        "Natural Jacobian requires a square cell block.");

    constexpr PetscInt blockSize =
        static_cast<PetscInt>(Indices::numPrimaryVariables);

    NaturalJacobianPattern pattern;
    pattern.blockSize = blockSize;

    const PetscInt ownedCells = backend.ownedCellCount(blockSize);
    pattern.rowBlocks.assign(static_cast<std::size_t>(ownedCells), PetscInt(-1));
    pattern.diagonalNnz.assign(static_cast<std::size_t>(ownedCells), 0);
    pattern.offDiagonalNnz.assign(static_cast<std::size_t>(ownedCells), 0);

    std::vector<std::vector<PetscInt>> rowColumns(
        static_cast<std::size_t>(ownedCells));

    std::unordered_map<PetscInt, std::vector<PetscInt>> extraCoupling;

    const auto map = backend.dofMap(blockSize);

    if constexpr (Indices::hasWellUnknown)
    {
        std::size_t perforationCount = 0;
        for (const auto &well : wells)
        {
            well.validate(backend.cellCount());
            perforationCount += well.perforations.size();
        }

        // StructuredGrid 的 DMDA global numbering 只能从 owned/ghost local map
        // 查询。井可能跨越任意 rank，因此 setup 时由各 owner 解析 block id，
        // 再一次 Allreduce 广播全部井拓扑；Newton 热循环不再做这类解析。
        std::vector<PetscInt> localWellBlocks(wells.size(), PetscInt(-1));
        std::vector<PetscInt> globalWellBlocks(wells.size(), PetscInt(-1));
        std::vector<PetscInt> localPerforationBlocks(perforationCount, PetscInt(-1));
        std::vector<PetscInt> globalPerforationBlocks(perforationCount, PetscInt(-1));

        std::size_t flat = 0;
        for (std::size_t wellIndex = 0; wellIndex < wells.size(); ++wellIndex)
        {
            const auto &well = wells[wellIndex];
            if (backend.isLocal(well.bhpCellId))
                localWellBlocks[wellIndex] = map.globalIndex(well.bhpCellId, 0) / blockSize;

            for (const auto &perforation : well.perforations)
            {
                if (backend.isLocal(perforation.currentCellId))
                    localPerforationBlocks[flat] =
                        map.globalIndex(perforation.currentCellId, 0) / blockSize;
                ++flat;
            }
        }

        if (!wells.empty())
        {
            PetscCallMPIAbort(
                backend.communicator(),
                MPI_Allreduce(
                    localWellBlocks.data(),
                    globalWellBlocks.data(),
                    static_cast<int>(localWellBlocks.size()),
                    MPIU_INT,
                    MPI_MAX,
                    backend.communicator()));
        }
        if (perforationCount > 0)
        {
            PetscCallMPIAbort(
                backend.communicator(),
                MPI_Allreduce(
                    localPerforationBlocks.data(),
                    globalPerforationBlocks.data(),
                    static_cast<int>(localPerforationBlocks.size()),
                    MPIU_INT,
                    MPI_MAX,
                    backend.communicator()));
        }

        flat = 0;
        for (std::size_t wellIndex = 0; wellIndex < wells.size(); ++wellIndex)
        {
            const auto &well = wells[wellIndex];
            const PetscInt bhpBlock = globalWellBlocks[wellIndex];
            if (bhpBlock < 0)
                throw std::runtime_error(
                    "Natural Jacobian pattern could not resolve a well BHP block.");

            if (backend.isLocal(well.bhpCellId))
            {
                auto &columns = extraCoupling[well.bhpCellId];
                for (std::size_t p = 0; p < well.perforations.size(); ++p)
                {
                    const PetscInt perforationBlock = globalPerforationBlocks[flat + p];
                    if (perforationBlock < 0)
                        throw std::runtime_error(
                            "Natural Jacobian pattern could not resolve a perforation block.");
                    columns.push_back(perforationBlock);
                }
            }

            for (const auto &perforation : well.perforations)
            {
                const PetscInt perforationBlock = globalPerforationBlocks[flat];
                if (perforationBlock < 0)
                    throw std::runtime_error(
                        "Natural Jacobian pattern could not resolve a perforation block.");
                if (backend.isLocal(perforation.currentCellId))
                    extraCoupling[perforation.currentCellId].push_back(bhpBlock);
                ++flat;
            }
        }
    }
    else if (!wells.empty())
    {
        throw std::invalid_argument(
            "Well specifications require Indices::hasWellUnknown.");
    }

    const PetscInt ownedDofBegin = backend.ownedDofBegin(blockSize);
    if (ownedDofBegin < 0 || ownedDofBegin % blockSize != 0)
        throw std::runtime_error(
            "Natural owned DOF range is not aligned to block rows.");

    const PetscInt firstOwnedBlock = ownedDofBegin / blockSize;
    const PetscInt endOwnedBlock = firstOwnedBlock + ownedCells;

    const auto globalBlock = [&](PetscInt cellId)
    {
        const PetscInt firstGlobal = map.globalIndex(cellId, 0);
        if (firstGlobal < 0 || firstGlobal % blockSize != 0)
            throw std::runtime_error(
                "Natural global DOF is not aligned to the Jacobian block size.");
        return firstGlobal / blockSize;
    };

    for (PetscInt cellId : backend.ownedCells())
    {
        const PetscInt localRow = backend.ownedBlockOrdinal(cellId, blockSize);
        if (localRow < 0 || localRow >= ownedCells)
            throw std::runtime_error(
                "Natural backend ownership is inconsistent with PETSc block rows.");

        auto &columns = rowColumns[static_cast<std::size_t>(localRow)];
        columns.push_back(globalBlock(cellId));

        for (PetscInt neighborId : backend.neighborCellIds(cellId))
            columns.push_back(globalBlock(neighborId));

        const auto extraIt = extraCoupling.find(cellId);
        if (extraIt != extraCoupling.end())
        {
            for (PetscInt columnBlock : extraIt->second)
                columns.push_back(columnBlock);
        }

        std::sort(columns.begin(), columns.end());
        columns.erase(std::unique(columns.begin(), columns.end()), columns.end());

        const PetscInt rowBlock = globalBlock(cellId);
        pattern.rowBlocks[static_cast<std::size_t>(localRow)] = rowBlock;

        PetscInt diagonal = 0;
        PetscInt offDiagonal = 0;
        for (PetscInt columnBlock : columns)
        {
            if (columnBlock >= firstOwnedBlock && columnBlock < endOwnedBlock)
                ++diagonal;
            else
                ++offDiagonal;
        }
        pattern.diagonalNnz[static_cast<std::size_t>(localRow)] = diagonal;
        pattern.offDiagonalNnz[static_cast<std::size_t>(localRow)] = offDiagonal;
    }

    pattern.rowOffsets.resize(static_cast<std::size_t>(ownedCells) + 1, 0);
    for (PetscInt localRow = 0; localRow < ownedCells; ++localRow)
    {
        if (pattern.rowBlocks[static_cast<std::size_t>(localRow)] < 0)
            throw std::runtime_error(
                "Natural Jacobian pattern is missing an owned block row.");

        pattern.rowOffsets[static_cast<std::size_t>(localRow + 1)] =
            pattern.rowOffsets[static_cast<std::size_t>(localRow)] +
            static_cast<PetscInt>(rowColumns[static_cast<std::size_t>(localRow)].size());
    }

    pattern.columnBlocks.reserve(
        static_cast<std::size_t>(pattern.rowOffsets.back()));
    for (const auto &columns : rowColumns)
        pattern.columnBlocks.insert(
            pattern.columnBlocks.end(), columns.begin(), columns.end());

    return pattern;
}

/**
 * @brief 按已缓存 pattern 创建固定稀疏结构的 Natural Jacobian。
 *
 * 标准 BAIJ 路径直接把精确 block-CSR 交给 PETSc 的 CSR preallocation API，
 * 由 PETSc 一次建立并组装全部 nonzero locations；只有用户通过 `-mat_type`
 * 选择其他 XAIJ 格式时才退回 G8J 的零 block 结构物化。随后固定 pattern，
 * Newton 阶段只更新既有数值 slot；MPI BAIJ 还会一次预热全部潜在远程 row
 * 目标，使后续 FINAL_ASSEMBLY 复用静态 off-process 通信图。
 */
template <class Backend>
[[nodiscard]] Mat createNaturalJacobianFromPattern(
    const Backend &backend,
    const NaturalJacobianPattern &pattern)
{
    if (pattern.blockSize <= 0)
        throw std::invalid_argument("Natural Jacobian pattern has invalid block size.");

    const PetscInt ownedCells = backend.ownedCellCount(pattern.blockSize);
    if (pattern.rowBlocks.size() != static_cast<std::size_t>(ownedCells) ||
        pattern.diagonalNnz.size() != static_cast<std::size_t>(ownedCells) ||
        pattern.offDiagonalNnz.size() != static_cast<std::size_t>(ownedCells) ||
        pattern.rowOffsets.size() != static_cast<std::size_t>(ownedCells) + 1)
        throw std::invalid_argument(
            "Natural Jacobian pattern does not match backend ownership.");

    Mat matrix = nullptr;
    const MPI_Comm comm = backend.communicator();

    PetscCallAbort(comm, MatCreate(comm, &matrix));
    PetscCallAbort(
        comm,
        MatSetSizes(
            matrix,
            backend.ownedDofCount(pattern.blockSize),
            backend.ownedDofCount(pattern.blockSize),
            backend.globalDofCount(pattern.blockSize),
            backend.globalDofCount(pattern.blockSize)));
    PetscCallAbort(comm, MatSetBlockSize(matrix, pattern.blockSize));
    PetscCallAbort(comm, MatSetType(matrix, MATBAIJ));
    PetscCallAbort(comm, MatSetFromOptions(matrix));

    PetscBool isSeqBaij = PETSC_FALSE;
    PetscBool isMpiBaij = PETSC_FALSE;
    PetscCallAbort(
        comm,
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(matrix), MATSEQBAIJ, &isSeqBaij));
    PetscCallAbort(
        comm,
        PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(matrix), MATMPIBAIJ, &isMpiBaij));

    /*
     * G8K: 对标准 BAIJ 路径直接把 G8J 已经构造好的精确 block-CSR 交给
     * PETSc。Mat{Seq,MPI}BAIJSetPreallocationCSR 不只分配容量，也一次建立
     * 精确 nonzero locations；因此不再需要逐 owned row 写零 dense block
     * 来“物化”结构。若用户通过 -mat_type 选择其他格式，仍保留通用 XAIJ
     * fallback，避免限制现有运行参数。
     */
    if (isSeqBaij)
    {
        PetscCallAbort(
            comm,
            MatSeqBAIJSetPreallocationCSR(
                matrix,
                pattern.blockSize,
                pattern.rowOffsets.data(),
                pattern.columnBlocks.data(),
                nullptr));
    }
    else if (isMpiBaij)
    {
        PetscCallAbort(
            comm,
            MatMPIBAIJSetPreallocationCSR(
                matrix,
                pattern.blockSize,
                pattern.rowOffsets.data(),
                pattern.columnBlocks.data(),
                nullptr));
    }
    else
    {
        PetscCallAbort(
            comm,
            MatXAIJSetPreallocation(
                matrix,
                pattern.blockSize,
                pattern.diagonalNnz.empty() ? nullptr : pattern.diagonalNnz.data(),
                pattern.offDiagonalNnz.empty() ? nullptr : pattern.offDiagonalNnz.data(),
                nullptr,
                nullptr));
        PetscCallAbort(
            comm,
            MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
        PetscCallAbort(comm, MatSetOption(matrix, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE));
        PetscCallAbort(comm, MatSetUp(matrix));

        std::size_t maximumColumnCount = 0;
        for (PetscInt localRow = 0; localRow < ownedCells; ++localRow)
        {
            const PetscInt begin = pattern.rowOffsets[static_cast<std::size_t>(localRow)];
            const PetscInt endOffset =
                pattern.rowOffsets[static_cast<std::size_t>(localRow + 1)];
            maximumColumnCount = std::max(
                maximumColumnCount,
                static_cast<std::size_t>(endOffset - begin));
        }

        const std::size_t blockValueCount =
            static_cast<std::size_t>(pattern.blockSize) *
            static_cast<std::size_t>(pattern.blockSize);
        std::vector<PetscScalar> zeroBlocks(
            maximumColumnCount * blockValueCount, PetscScalar(0.0));

        for (PetscInt localRow = 0; localRow < ownedCells; ++localRow)
        {
            const PetscInt begin = pattern.rowOffsets[static_cast<std::size_t>(localRow)];
            const PetscInt endOffset =
                pattern.rowOffsets[static_cast<std::size_t>(localRow + 1)];
            const PetscInt columnCount = endOffset - begin;
            if (columnCount <= 0)
                continue;

            const PetscInt rowBlock =
                pattern.rowBlocks[static_cast<std::size_t>(localRow)];
            PetscCallAbort(
                comm,
                MatSetValuesBlocked(
                    matrix,
                    1,
                    &rowBlock,
                    columnCount,
                    pattern.columnBlocks.data() + begin,
                    zeroBlocks.data(),
                    INSERT_VALUES));
        }

        PetscCallAbort(comm, MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY));
        PetscCallAbort(comm, MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY));
    }

    PetscCallAbort(comm, MatSetOption(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
    PetscCallAbort(comm, MatSetOption(matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));

    /*
     * Column-oriented AD 会把 reciprocal row 写到邻居 owner，因此 Newton 的
     * FINAL_ASSEMBLY 仍是必要的：PETSc 必须把 off-process stash 发送给真实
     * row owner。不能安全删除这一步，也不能设置 MAT_NO_OFF_PROC_ENTRIES。
     *
     * 但 off-process 目标完全由静态 sparsity 决定。这里在 setup 时把每个
     * local-row/offdiag-column pair 的反向 block (remote row, local column)
     * 以零值发送一次；它覆盖 TPFA reciprocal 以及 well BHP/RATE 的全部潜在
     * 远程目标。随后 MAT_SUBSET_OFF_PROC_ENTRIES 可让后续 assembly 复用固定
     * 邻居通信图，跳过每次 rendezvous。该代价只支付一次，不进入 Newton。
     */
    int communicatorSize = 1;
    PetscCallMPIAbort(comm, MPI_Comm_size(comm, &communicatorSize));
    if (communicatorSize > 1)
    {
        const PetscInt ownedDofBegin = backend.ownedDofBegin(pattern.blockSize);
        if (ownedDofBegin < 0 || ownedDofBegin % pattern.blockSize != 0)
            throw std::runtime_error(
                "Natural Jacobian ownership is not aligned to block rows.");

        const PetscInt firstOwnedBlock = ownedDofBegin / pattern.blockSize;
        const PetscInt endOwnedBlock = firstOwnedBlock + ownedCells;
        const std::size_t blockValueCount =
            static_cast<std::size_t>(pattern.blockSize) *
            static_cast<std::size_t>(pattern.blockSize);
        std::vector<PetscScalar> zeroBlock(blockValueCount, PetscScalar(0.0));

        PetscCallAbort(
            comm,
            MatSetOption(matrix, MAT_SUBSET_OFF_PROC_ENTRIES, PETSC_TRUE));

        for (PetscInt localRow = 0; localRow < ownedCells; ++localRow)
        {
            const PetscInt localColumnBlock =
                pattern.rowBlocks[static_cast<std::size_t>(localRow)];
            const PetscInt begin = pattern.rowOffsets[static_cast<std::size_t>(localRow)];
            const PetscInt endOffset =
                pattern.rowOffsets[static_cast<std::size_t>(localRow + 1)];

            for (PetscInt offset = begin; offset < endOffset; ++offset)
            {
                const PetscInt remoteRowBlock =
                    pattern.columnBlocks[static_cast<std::size_t>(offset)];
                if (remoteRowBlock >= firstOwnedBlock && remoteRowBlock < endOwnedBlock)
                    continue;

                PetscCallAbort(
                    comm,
                    MatSetValuesBlocked(
                        matrix,
                        1,
                        &remoteRowBlock,
                        1,
                        &localColumnBlock,
                        zeroBlock.data(),
                        ADD_VALUES));
            }
        }

        PetscCallAbort(comm, MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY));
        PetscCallAbort(comm, MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY));
    }

    return matrix;
}

/** @brief 兼容独立调用：构造 pattern 后创建矩阵。 */
template <class Indices, class Backend>
[[nodiscard]] Mat createNaturalJacobian(
    const Backend &backend,
    const std::vector<NaturalWell<Indices>> &wells)
{
    const auto pattern = buildNaturalJacobianPattern<Indices>(backend, wells);
    return createNaturalJacobianFromPattern(backend, pattern);
}

} // namespace MPMC
