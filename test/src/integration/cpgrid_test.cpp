/**
 * @file cpgrid_test.cpp
 * @brief 集成测试：验证 `cpgrid_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <cpgrid/cpgrid.hpp>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

/**
 * @brief 本测试只要求主状态每单元 4 个自由度。
 */
struct ValidationTag
{
    static constexpr PetscInt numVars_ = 4;
};

constexpr PetscInt kDof =
    ValidationTag::numVars_;

constexpr PetscReal kTolerance =
    1.0e-12;

/**
 * @brief 测试结果。
 */
struct TestResult
{
    std::string name;
    PetscReal maxError{0.0};
    bool passed{false};
};

/**
 * @brief 对所有 MPI rank 的最大误差取全局最大值。
 */
PetscReal globalMax(
    MPI_Comm comm,
    PetscReal localValue)
{
    PetscReal globalValue = 0.0;

    PetscCallMPIAbort(
        comm,
        MPI_Allreduce(
            &localValue,
            &globalValue,
            1,
            MPIU_REAL,
            MPI_MAX,
            comm));

    return globalValue;
}

/**
 * @brief 为 global Vec 填入一个由全局 DOF 唯一确定的解析值。
 *
 * 使用 0.25 + globalDof，便于同时检查 owned 和 ghost 数据。
 */
void fillGlobalReference(
    Vec global)
{
    MPI_Comm comm =
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                global));

    PetscInt begin = 0;
    PetscInt end = 0;

    PetscCallAbort(
        comm,
        VecGetOwnershipRange(
            global,
            &begin,
            &end));

    PetscScalar *array = nullptr;

    PetscCallAbort(
        comm,
        VecGetArray(
            global,
            &array));

    for (PetscInt globalDof = begin;
         globalDof < end;
         ++globalDof)
    {
        array[globalDof - begin] =
            static_cast<PetscScalar>(
                0.25 +
                static_cast<PetscReal>(
                    globalDof));
    }

    PetscCallAbort(
        comm,
        VecRestoreArray(
            global,
            &array));
}

/**
 * @brief 验证 global -> local scatter。
 *
 * 本地快照的每个位置都通过 DofMap::localToGlobalIndices() 得到其期望全局 DOF，
 * 因而可以直接检查 ghost 数据是否来自正确的远端位置。
 */
TestResult testGlobalToLocal(
    const MPMC::CpGrid<ValidationTag> &grid)
{
    const auto &map =
        grid.dofMap(kDof);

    Vec global =
        grid.createGlobalVector(kDof);

    fillGlobalReference(global);

    Vec local =
        grid.createLocalVector(
            kDof,
            global);

    PetscInt localSize = 0;

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecGetLocalSize(
            local,
            &localSize));

    if (localSize !=
        map.localSnapshotDofCount())
    {
        PetscCallAbort(
            PETSC_COMM_SELF,
            VecDestroy(&local));

        PetscCallAbort(
            grid.mesh().communicator(),
            VecDestroy(&global));

        throw std::runtime_error(
            "global->local test: local snapshot size is inconsistent with DofMap.");
    }

    const auto localToGlobal =
        map.localToGlobalIndices();

    if (static_cast<PetscInt>(
            localToGlobal.size()) !=
        localSize)
    {
        PetscCallAbort(
            PETSC_COMM_SELF,
            VecDestroy(&local));

        PetscCallAbort(
            grid.mesh().communicator(),
            VecDestroy(&global));

        throw std::runtime_error(
            "global->local test: local-to-global index count is inconsistent.");
    }

    const PetscScalar *array = nullptr;

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecGetArrayRead(
            local,
            &array));

    PetscReal localMaxError = 0.0;

    for (PetscInt localDof = 0;
         localDof < localSize;
         ++localDof)
    {
        const PetscInt globalDof =
            localToGlobal[
                static_cast<std::size_t>(
                    localDof)];

        const PetscScalar expected =
            static_cast<PetscScalar>(
                0.25 +
                static_cast<PetscReal>(
                    globalDof));

        localMaxError =
            std::max(
                localMaxError,
                PetscAbsScalar(
                    array[localDof] -
                    expected));
    }

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecRestoreArrayRead(
            local,
            &array));

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecDestroy(&local));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecDestroy(&global));

    const PetscReal error =
        globalMax(
            grid.mesh().communicator(),
            localMaxError);

    return {
        "global -> local owned+ghost scatter",
        error,
        error <= kTolerance};
}

/**
 * @brief 验证 local -> global 的 owned 数据回写。
 *
 * 仅给 owned 区写入非零值，ghost 区保持 0，因此每个 global DOF 应只收到
 * owner rank 的唯一贡献。
 */
TestResult testLocalToGlobalOwned(
    const MPMC::CpGrid<ValidationTag> &grid)
{
    const auto &map =
        grid.dofMap(kDof);

    Vec local =
        grid.createLocalVector(kDof);

    Vec global =
        grid.createGlobalVector(kDof);

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecSet(
            local,
            static_cast<PetscScalar>(0.0)));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecSet(
            global,
            static_cast<PetscScalar>(0.0)));

    PetscScalar *localArray = nullptr;

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecGetArray(
            local,
            &localArray));

    for (const auto &cell :
         grid.localCells())
    {
        for (PetscInt component = 0;
             component < kDof;
             ++component)
        {
            const PetscInt localDof =
                map.localIndex(
                    cell,
                    component);

            localArray[localDof] =
                static_cast<PetscScalar>(
                    10.0 +
                    static_cast<PetscReal>(
                        component));
        }
    }

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecRestoreArray(
            local,
            &localArray));

    grid.localToGlobal(
        kDof,
        local,
        global,
        ADD_VALUES);

    PetscInt begin = 0;
    PetscInt end = 0;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetOwnershipRange(
            global,
            &begin,
            &end));

    if (begin !=
            map.firstOwnedGlobalDof() ||
        end !=
            map.endOwnedGlobalDof())
    {
        PetscCallAbort(
            PETSC_COMM_SELF,
            VecDestroy(&local));

        PetscCallAbort(
            grid.mesh().communicator(),
            VecDestroy(&global));

        throw std::runtime_error(
            "local->global owned test: PETSc ownership range differs from DofMap.");
    }

    const PetscScalar *globalArray =
        nullptr;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetArrayRead(
            global,
            &globalArray));

    PetscReal localMaxError = 0.0;

    for (const auto &cell :
         grid.localCells())
    {
        for (PetscInt component = 0;
             component < kDof;
             ++component)
        {
            const PetscInt globalDof =
                map.globalIndex(
                    cell,
                    component);

            const PetscScalar expected =
                static_cast<PetscScalar>(
                    10.0 +
                    static_cast<PetscReal>(
                        component));

            const PetscScalar actual =
                globalArray[
                    globalDof -
                    begin];

            localMaxError =
                std::max(
                    localMaxError,
                    PetscAbsScalar(
                        actual -
                        expected));
        }
    }

    PetscCallAbort(
        grid.mesh().communicator(),
        VecRestoreArrayRead(
            global,
            &globalArray));

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecDestroy(&local));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecDestroy(&global));

    const PetscReal error =
        globalMax(
            grid.mesh().communicator(),
            localMaxError);

    return {
        "local -> global owned contribution",
        error,
        error <= kTolerance};
}

/**
 * @brief 收集所有 rank 的 ghost current-cell ids，并计算每个 cell 被多少个远端
 * rank 作为 ghost 使用。
 */
std::vector<PetscInt>
computeGhostMultiplicity(
    const MPMC::DofMap &map)
{
    const MPMC::Mesh &mesh =
        map.mesh();

    const auto &localGhosts =
        map.ghostCellIds();

    if (localGhosts.size() >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max()))
    {
        throw std::overflow_error(
            "Local ghost count exceeds MPI_Allgatherv int count.");
    }

    const int localCount =
        static_cast<int>(
            localGhosts.size());

    std::vector<int> counts(
        static_cast<std::size_t>(
            mesh.processCount()),
        0);

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Allgather(
            &localCount,
            1,
            MPI_INT,
            counts.data(),
            1,
            MPI_INT,
            mesh.communicator()));

    std::vector<int> displacements(
        counts.size(),
        0);

    int totalCount = 0;

    for (std::size_t rank = 0;
         rank < counts.size();
         ++rank)
    {
        displacements[rank] =
            totalCount;

        if (counts[rank] >
            std::numeric_limits<int>::max() -
                totalCount)
        {
            throw std::overflow_error(
                "Global ghost count exceeds MPI_Allgatherv int count.");
        }

        totalCount +=
            counts[rank];
    }

    std::vector<PetscInt> allGhosts(
        static_cast<std::size_t>(
            totalCount));

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Allgatherv(
            localGhosts.empty()
                ? nullptr
                : localGhosts.data(),
            localCount,
            MPIU_INT,
            allGhosts.empty()
                ? nullptr
                : allGhosts.data(),
            counts.data(),
            displacements.data(),
            MPIU_INT,
            mesh.communicator()));

    std::vector<PetscInt> multiplicity(
        mesh.cellCount(),
        0);

    for (PetscInt cellId :
         allGhosts)
    {
        if (cellId < 0 ||
            static_cast<std::size_t>(
                cellId) >=
                multiplicity.size())
        {
            throw std::runtime_error(
                "Ghost current-cell id is outside Mesh range.");
        }

        ++multiplicity[
            static_cast<std::size_t>(
                cellId)];
    }

    return multiplicity;
}

/**
 * @brief 验证 ghost 区通过 local -> global ADD_VALUES 回写到真正 owner。
 *
 * 每个 rank 把自己的所有 ghost cell block 置为 1，owned block 置为 0。
 * 因此某个 global cell 最终应等于：
 *
 *   “该 cell 在多少个远端 rank 上作为 ghost 出现”
 */
TestResult testLocalToGlobalGhost(
    const MPMC::CpGrid<ValidationTag> &grid)
{
    const auto &map =
        grid.dofMap(kDof);

    const auto multiplicity =
        computeGhostMultiplicity(map);

    Vec local =
        grid.createLocalVector(kDof);

    Vec global =
        grid.createGlobalVector(kDof);

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecSet(
            local,
            static_cast<PetscScalar>(0.0)));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecSet(
            global,
            static_cast<PetscScalar>(0.0)));

    PetscScalar *localArray = nullptr;

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecGetArray(
            local,
            &localArray));

    for (PetscInt ghostCellId :
         map.ghostCellIds())
    {
        const auto &ghostCell =
            grid.mesh().cellByCurrentId(
                ghostCellId);

        for (PetscInt component = 0;
             component < kDof;
             ++component)
        {
            localArray[
                map.localIndex(
                    ghostCell,
                    component)] =
                static_cast<PetscScalar>(
                    1.0);
        }
    }

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecRestoreArray(
            local,
            &localArray));

    grid.localToGlobal(
        kDof,
        local,
        global,
        ADD_VALUES);

    PetscInt begin = 0;
    PetscInt end = 0;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetOwnershipRange(
            global,
            &begin,
            &end));

    const PetscScalar *globalArray =
        nullptr;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetArrayRead(
            global,
            &globalArray));

    PetscReal localMaxError = 0.0;

    for (const auto &cell :
         grid.localCells())
    {
        const PetscScalar expected =
            static_cast<PetscScalar>(
                multiplicity[
                    static_cast<std::size_t>(
                        cell.id())]);

        for (PetscInt component = 0;
             component < kDof;
             ++component)
        {
            const PetscInt globalDof =
                map.globalIndex(
                    cell,
                    component);

            const PetscScalar actual =
                globalArray[
                    globalDof -
                    begin];

            localMaxError =
                std::max(
                    localMaxError,
                    PetscAbsScalar(
                        actual -
                        expected));
        }
    }

    PetscCallAbort(
        grid.mesh().communicator(),
        VecRestoreArrayRead(
            global,
            &globalArray));

    PetscCallAbort(
        PETSC_COMM_SELF,
        VecDestroy(&local));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecDestroy(&global));

    const PetscReal error =
        globalMax(
            grid.mesh().communicator(),
            localMaxError);

    return {
        "local -> global ghost ADD_VALUES",
        error,
        error <= kTolerance};
}

/**
 * @brief 构造一个 bs×bs 的对角块。
 */
std::array<PetscScalar, kDof * kDof>
diagonalBlock(
    PetscScalar diagonalValue)
{
    std::array<
        PetscScalar,
        kDof * kDof>
        block{};

    for (PetscInt component = 0;
         component < kDof;
         ++component)
    {
        block[
            static_cast<std::size_t>(
                component * kDof +
                component)] =
            diagonalValue;
    }

    return block;
}

/**
 * @brief 验证 blocked-local Jacobian 装配、跨 rank 列映射和矩阵乘法。
 *
 * 每个 owned cell row 装配：
 *
 *   A_ii = 10 I
 *   A_ij = -I   对每个内部面邻居 j
 *
 * 然后令 x=1，则每个分量解析结果为：
 *
 *   y_i = 10 - N_internal_faces(i)
 */
TestResult testBlockedMatrix(
    const MPMC::CpGrid<ValidationTag> &grid,
    PetscLogDouble &matrixMallocs)
{
    const auto &map =
        grid.dofMap(kDof);

    Mat matrix =
        grid.createMatrix(kDof);

    const auto diagonal =
        diagonalBlock(
            static_cast<PetscScalar>(
                10.0));

    const auto offDiagonal =
        diagonalBlock(
            static_cast<PetscScalar>(
                -1.0));

    for (const auto &cell :
         grid.localCells())
    {
        const PetscInt localRow =
            map.localBlockIndex(cell);

        PetscCallAbort(
            grid.mesh().communicator(),
            MatSetValuesBlockedLocal(
                matrix,
                1,
                &localRow,
                1,
                &localRow,
                diagonal.data(),
                ADD_VALUES));

        for (const auto &face :
             cell.faces())
        {
            const auto *neighbor =
                face.neighborCell();

            if (neighbor == nullptr)
            {
                continue;
            }

            const PetscInt localColumn =
                map.localBlockIndex(
                    *neighbor);

            PetscCallAbort(
                grid.mesh().communicator(),
                MatSetValuesBlockedLocal(
                    matrix,
                    1,
                    &localRow,
                    1,
                    &localColumn,
                    offDiagonal.data(),
                    ADD_VALUES));
        }
    }

    PetscCallAbort(
        grid.mesh().communicator(),
        MatAssemblyBegin(
            matrix,
            MAT_FINAL_ASSEMBLY));

    PetscCallAbort(
        grid.mesh().communicator(),
        MatAssemblyEnd(
            matrix,
            MAT_FINAL_ASSEMBLY));

    MatInfo info{};

    PetscCallAbort(
        grid.mesh().communicator(),
        MatGetInfo(
            matrix,
            MAT_LOCAL,
            &info));

    matrixMallocs =
        info.mallocs;

    Vec x =
        grid.createGlobalVector(kDof);

    Vec y =
        grid.createGlobalVector(kDof);

    PetscCallAbort(
        grid.mesh().communicator(),
        VecSet(
            x,
            static_cast<PetscScalar>(
                1.0)));

    PetscCallAbort(
        grid.mesh().communicator(),
        MatMult(
            matrix,
            x,
            y));

    const PetscScalar *yArray =
        nullptr;

    PetscCallAbort(
        grid.mesh().communicator(),
        VecGetArrayRead(
            y,
            &yArray));

    PetscReal localMaxError = 0.0;

    for (const auto &cell :
         grid.localCells())
    {
        PetscInt internalFaceCount = 0;

        for (const auto &face :
             cell.faces())
        {
            if (!face.isBoundary())
            {
                ++internalFaceCount;
            }
        }

        const PetscScalar expected =
            static_cast<PetscScalar>(
                10.0 -
                static_cast<PetscReal>(
                    internalFaceCount));

        const PetscInt localBlock =
            map.localBlockIndex(cell);

        for (PetscInt component = 0;
             component < kDof;
             ++component)
        {
            const PetscInt localDof =
                localBlock *
                    kDof +
                component;

            localMaxError =
                std::max(
                    localMaxError,
                    PetscAbsScalar(
                        yArray[localDof] -
                        expected));
        }
    }

    PetscCallAbort(
        grid.mesh().communicator(),
        VecRestoreArrayRead(
            y,
            &yArray));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecDestroy(&x));

    PetscCallAbort(
        grid.mesh().communicator(),
        VecDestroy(&y));

    PetscCallAbort(
        grid.mesh().communicator(),
        MatDestroy(&matrix));

    const PetscReal error =
        globalMax(
            grid.mesh().communicator(),
            localMaxError);

    return {
        "MatSetValuesBlockedLocal + assembly + MatMult",
        error,
        error <= kTolerance};
}

/**
 * @brief 打印每个 rank 的布局摘要。
 */
void printLayoutSummary(
    const MPMC::CpGrid<ValidationTag> &grid)
{
    const auto &map =
        grid.dofMap(kDof);

    PetscSynchronizedPrintf(
        grid.mesh().communicator(),
        "[rank %d] ownedCells=%"
        PetscInt_FMT
        ", ghostCells=%"
        PetscInt_FMT
        ", ownedDOFs=%"
        PetscInt_FMT
        ", localSnapshotDOFs=%"
        PetscInt_FMT
        ", firstOwnedCell=%"
        PetscInt_FMT
        "\n",
        grid.mesh().rank(),
        map.ownedCellCount(),
        map.ghostCellCount(),
        map.ownedDofCount(),
        map.localSnapshotDofCount(),
        map.firstOwnedCellId());

    PetscSynchronizedFlush(
        grid.mesh().communicator(),
        PETSC_STDOUT);
}

/**
 * @brief 打印单项测试结果。
 */
void printResult(
    MPI_Comm comm,
    const TestResult &result)
{
    PetscPrintf(
        comm,
        "[%s] %-48s maxError = %.3e\n",
        result.passed
            ? "PASS"
            : "FAIL",
        result.name.c_str(),
        static_cast<double>(
            result.maxError));
}

void testGrdeclToleranceNodeMerging()
{
    MPMC::GrdeclGridData data;
    data.nx = 2;
    data.ny = 1;
    data.nz = 1;
    data.cartesianCellCount = 2;
    data.nodeMergeTolerance = 1.0e-3;
    data.actnum = {1, 1};
    data.activeCells.resize(2);

    const auto fillCell = [](
                              MPMC::GrdeclCell &cell,
                              int cartesianId,
                              double xMinus,
                              double xPlus) {
        cell.cartesianId = cartesianId;
        cell.i = cartesianId;
        for (int iz = 0; iz < 2; ++iz)
            for (int iy = 0; iy < 2; ++iy)
                for (int ix = 0; ix < 2; ++ix)
                {
                    const auto corner = static_cast<std::size_t>(
                        ix + 2 * iy + 4 * iz);
                    cell.corner[corner] = {
                        ix == 0 ? xMinus : xPlus,
                        static_cast<double>(iy),
                        static_cast<double>(iz)};
                }
    };

    // The two shared-face x coordinates straddle an llround() bucket
    // boundary but are much closer than the requested merge tolerance.
    fillCell(data.activeCells[0], 0, 0.0, 1.0004999);
    fillCell(data.activeCells[1], 1, 1.0005001, 2.0);
    data.activeCells[0].neighborStorage[1] = 1;
    data.activeCells[1].neighborStorage[0] = 0;

    MPMC::Mesh mesh(data);
    if (mesh.nodeCount() != 12)
        throw std::runtime_error(
            "GRDECL tolerance merge must share the four conforming face nodes.");

    MPMC::GrdeclGridData diagonal = data;
    diagonal.activeCells = {};
    diagonal.activeCells.resize(2);
    fillCell(diagonal.activeCells[0], 0, 0.0, 1.0);
    fillCell(diagonal.activeCells[1], 1, 3.0, 4.0);
    for (auto &corner : diagonal.activeCells[1].corner)
    {
        corner.y += 3.0;
        corner.z += 3.0;
    }
    diagonal.activeCells[0].corner[0] = {-0.00049, -0.00049, -0.00049};
    diagonal.activeCells[1].corner[0] = {0.00049, 0.00049, 0.00049};

    MPMC::Mesh diagonalMesh(diagonal);
    if (diagonalMesh.nodeCount() != 16)
        throw std::runtime_error(
            "GRDECL node merge must use true distance, not quantized-key equality.");

    MPMC::GrdeclGridData flow = data;
    flow.nodeMergeTolerance = 1.0e-6;
    flow.activeCells = {};
    flow.activeCells.resize(2);
    fillCell(flow.activeCells[0], 0, 0.0, 1.0);
    fillCell(flow.activeCells[1], 1, 1.1, 2.1);
    for (auto &corner : flow.activeCells[1].corner)
    {
        corner.y *= 2.0;
        corner.z *= 2.0;
    }
    flow.activeCells[0].neighborStorage[1] = 1;
    flow.activeCells[1].neighborStorage[0] = 0;

    MPMC::Mesh flowMesh(flow, PETSC_COMM_SELF);
    flowMesh.prepareForUse();
    MPMC::CpGrid<ValidationTag> flowGrid(flowMesh);
    PetscCallAbort(
        flowMesh.communicator(),
        VecSet(flowGrid.permeabilityVector(), 1.0));
    PetscCallAbort(
        flowMesh.communicator(),
        VecSet(flowGrid.porosityVector(), 0.2));
    flowGrid.setup();

    const auto &firstCell = flowMesh.cellByInputIndex(0);
    const auto face = std::find_if(
        firstCell.faces().begin(),
        firstCell.faces().end(),
        [](const MPMC::Face &candidate) {
            return candidate.inputPosition() == 2;
        });
    if (face == firstCell.faces().end() || face->neighborFace() == nullptr)
        throw std::runtime_error(
            "Internal GRDECL face must expose its opposite face instance.");

    // Half transmissibilities are 1/(0.5/1)=2 and 4/(0.5/1)=8;
    // the harmonic connection is therefore 1/(1/2+1/8)=1.6.
    if (std::abs(flowGrid.transmissibility(*face) - 1.6) > 1.0e-12)
        throw std::runtime_error(
            "CpGrid must use each cell's own half-face geometry.");
}

/**
 * @brief 执行全部并行基础设施测试。
 */
bool runValidation(
    const std::string &meshDirectory)
{
    testGrdeclToleranceNodeMerging();

    MPMC::Mesh mesh(
        meshDirectory);

    mesh.prepareForUse();

    MPMC::CpGrid<ValidationTag>
        grid(mesh);

    printLayoutSummary(grid);

    PetscInt totalGhostCells = 0;
    const PetscInt localGhostCells =
        grid.dofMap(kDof)
            .ghostCellCount();

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Allreduce(
            &localGhostCells,
            &totalGhostCells,
            1,
            MPIU_INT,
            MPI_SUM,
            mesh.communicator()));

    PetscPrintf(
        mesh.communicator(),
        "\nCpGrid parallel validation\n"
        "  ranks             = %d\n"
        "  cells             = %zu\n"
        "  dof/cell          = %"
        PetscInt_FMT
        "\n"
        "  global DOFs       = %"
        PetscInt_FMT
        "\n"
        "  total ghost cells = %"
        PetscInt_FMT
        "\n\n",
        mesh.processCount(),
        mesh.cellCount(),
        kDof,
        grid.dofMap(kDof)
            .globalDofCount(),
        totalGhostCells);

    const TestResult globalToLocal =
        testGlobalToLocal(grid);

    const TestResult ownedToGlobal =
        testLocalToGlobalOwned(grid);

    const TestResult ghostToGlobal =
        testLocalToGlobalGhost(grid);

    PetscLogDouble matrixMallocs = 0.0;

    const TestResult matrix =
        testBlockedMatrix(
            grid,
            matrixMallocs);

    printResult(
        mesh.communicator(),
        globalToLocal);

    printResult(
        mesh.communicator(),
        ownedToGlobal);

    printResult(
        mesh.communicator(),
        ghostToGlobal);

    printResult(
        mesh.communicator(),
        matrix);

    PetscPrintf(
        mesh.communicator(),
        "[INFO] matrix MatInfo.mallocs = %.0f\n",
        static_cast<double>(
            matrixMallocs));

    if (mesh.processCount() > 1 &&
        totalGhostCells == 0)
    {
        PetscPrintf(
            mesh.communicator(),
            "[WARN] multi-rank run produced zero ghost cells; "
            "the mesh may be disconnected or partitioning is unusual.\n");
    }

    const bool passed =
        globalToLocal.passed &&
        ownedToGlobal.passed &&
        ghostToGlobal.passed &&
        matrix.passed;

    PetscPrintf(
        mesh.communicator(),
        "\n============================================\n"
        "CpGrid parallel validation: %s\n"
        "============================================\n",
        passed
            ? "ALL PASS"
            : "FAILED");

    return passed;
}

} // namespace

int main(
    int argc,
    char **argv)
{
    PetscInitialize(
        &argc,
        &argv,
        nullptr,
        nullptr);

    int exitCode = 0;

    try
    {
        char meshDirectory[
            PETSC_MAX_PATH_LEN] = {};

        PetscBool found =
            PETSC_FALSE;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetString(
                nullptr,
                nullptr,
                "-mesh_dir",
                meshDirectory,
                sizeof(meshDirectory),
                &found));

        if (!found)
        {
            throw std::invalid_argument(
                "Please provide -mesh_dir <MRST-export-directory>.");
        }

        /*
         * 所有 PETSc owner 在 PetscFinalize() 前析构。
         */
        const bool passed =
            runValidation(
                std::string(
                    meshDirectory));

        if (!passed)
        {
            exitCode = 2;
        }
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            PETSC_COMM_WORLD,
            "CpGrid validation failed with exception: %s\n",
            error.what());

        exitCode = 1;
    }

    PetscFinalize();
    return exitCode;
}
