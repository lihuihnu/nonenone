/**
 * @file structuredgrid_backend.hpp
 * @brief Natural 运行时面向 StructuredGrid 的网格后端。
 */
#pragma once

#include <natural/petsc/grid_backend_common.hpp>
#include <structuredgrid/structuredgrid.hpp>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscis.h>
#include <petscsys.h>
#include <petscvec.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace MPMC
{

/**
 * @brief StructuredGrid/DMDA 到 Natural PETSc runtime 的后端适配器。
 *
 * StructuredGrid 保留 Cartesian cell id；PETSc 全局 DOF 通过 DMDA 的
 * local-to-global mapping 获取，避免假设多进程下 global = cell*dof。
 */
template <class Indices, class Grid>
class NaturalStructuredGridBackend final
{
public:
    using CellId = PetscInt;
    using GridType = Grid;
    using DofMap = NaturalDofMapView<NaturalStructuredGridBackend>;

private:
    struct LayoutInfo final
    {
        PetscInt dof{0};
        DMDARegion owned{};
        DMDARegion ghost{};
        DMDAStencilType stencilType{DMDA_STENCIL_STAR};
        PetscInt ownedDofBegin{0};
        PetscInt ownedDofEnd{0};
        PetscInt globalDofCount{0};
        ISLocalToGlobalMapping localToGlobal{nullptr}; // borrowed from DM
    };

public:
    class LocalReadView final
    {
    public:
        LocalReadView(
            const NaturalStructuredGridBackend &backend,
            PetscInt dofPerCell,
            Vec global)
            : backend_(backend),
              dofPerCell_(dofPerCell),
              local_(backend_.grid_.borrowLocalVector(
                  static_cast<int>(dofPerCell_),
                  global))
        {
            PetscCallAbort(
                backend_.communicator(),
                VecGetArrayRead(local_, &array_));
        }

        LocalReadView(const LocalReadView &) = delete;
        LocalReadView &operator=(const LocalReadView &) = delete;
        LocalReadView(LocalReadView &&) = delete;
        LocalReadView &operator=(LocalReadView &&) = delete;

        ~LocalReadView() noexcept
        {
            if (local_ != nullptr && array_ != nullptr)
            {
                PetscCallAbort(
                    backend_.communicator(),
                    VecRestoreArrayRead(local_, &array_));
            }
            if (local_ != nullptr)
            {
                backend_.grid_.restoreLocalVector(
                    static_cast<int>(dofPerCell_),
                    local_);
            }
        }

        [[nodiscard]] double value(
            CellId cellId,
            PetscInt component) const
        {
            const PetscInt index =
                backend_.localIndex(
                    cellId,
                    dofPerCell_,
                    component);
            return PetscRealPart(array_[index]);
        }

        /** @brief 调用方已知局部 block 序号时快速访问局部向量。 */
        [[nodiscard]] double valueByLocalBlock(
            PetscInt localBlock,
            PetscInt component) const
        {
            if (localBlock < 0 || component < 0 || component >= dofPerCell_)
                throw std::out_of_range("Natural local Vec block/component is out of range.");
            return PetscRealPart(
                array_[localBlock * dofPerCell_ + component]);
        }

    private:
        const NaturalStructuredGridBackend &backend_;
        PetscInt dofPerCell_{0};
        Vec local_{nullptr};
        const PetscScalar *array_{nullptr};
    };

    class AssemblyView final
    {
    public:
        explicit AssemblyView(
            const NaturalStructuredGridBackend &backend)
            : backend_(backend),
              gridAccess_(backend_.grid_)
        {
        }

        AssemblyView(const AssemblyView &) = delete;
        AssemblyView &operator=(const AssemblyView &) = delete;
        AssemblyView(AssemblyView &&) = delete;
        AssemblyView &operator=(AssemblyView &&) = delete;

        ~AssemblyView() noexcept = default;

        [[nodiscard]] double porosity(CellId cellId) const
        {
            return backend_.grid_.porosity(
                static_cast<int>(cellId));
        }

        [[nodiscard]] double cellVolume(CellId cellId) const
        {
            return backend_.grid_.cellVolume(
                static_cast<int>(cellId));
        }

        [[nodiscard]] std::vector<NaturalCellConnection>
        connections(CellId cellId) const
        {
            const auto faces = backend_.grid_.faces(
                static_cast<int>(cellId));

            std::vector<NaturalCellConnection> result;
            result.reserve(faces.size());

            for (const auto &face : faces)
            {
                const int neighbor = face.neighbor();

                result.push_back({
                    static_cast<PetscInt>(neighbor),
                    backend_.grid_.transmissibility(
                        static_cast<int>(cellId),
                        face),
                    backend_.grid_.gravityPotentialDifference(
                        static_cast<int>(cellId),
                        face)});
            }

            return result;
        }

    private:
        const NaturalStructuredGridBackend &backend_;
        typename Grid::AssemblyAccess gridAccess_;
    };

    explicit NaturalStructuredGridBackend(Grid &grid)
        : grid_(grid)
    {
        validateAllCellsActive_();
    }

    [[nodiscard]] Grid &grid() noexcept { return grid_; }
    [[nodiscard]] const Grid &grid() const noexcept { return grid_; }

    void registerLayout(PetscInt dof)
    {
        grid_.registerLayout(static_cast<int>(dof));
        layoutCache_.erase(dof);
    }

    /** @brief 保留历史 StructuredGrid AD 辅助布局，由 Natural 显式请求。 */
    void registerAdAuxiliaryLayout(PetscInt dof)
    {
        registerLayout(dof);
    }

    [[nodiscard]] Vec createGlobalVector(PetscInt dof) const
    {
        return grid_.createGlobalVector(static_cast<int>(dof));
    }

    /**
     * @brief 将求解器 Vec 复制为标准 Cartesian 单元顺序用于输出。
     *
     * DMDA may use a rank-dependent PETSc global numbering.  External files should
     * instead remain stable as `cell_id -> component`, where `cell_id` is the
     * StructuredGrid Cartesian global cell id.
     */
    [[nodiscard]] Vec createInputOrderedCopy(
        PetscInt dofPerCell,
        Vec currentOrderedVector) const
    {
        if (currentOrderedVector == nullptr)
            throw std::invalid_argument(
                "createInputOrderedCopy requires a valid PETSc Vec.");

        PetscInt globalSize = 0;
        PetscCallAbort(
            communicator(),
            VecGetSize(currentOrderedVector, &globalSize));

        if (globalSize != cellCount() * dofPerCell)
            throw std::invalid_argument(
                "Vector size does not match StructuredGrid cell layout.");

        PetscInt localSize = 0;
        PetscCallAbort(
            communicator(),
            VecGetLocalSize(currentOrderedVector, &localSize));

        Vec ordered = nullptr;
        PetscCallAbort(
            communicator(),
            VecCreateMPI(
                communicator(),
                localSize,
                globalSize,
                &ordered));
        PetscCallAbort(
            communicator(),
            VecSetBlockSize(ordered, dofPerCell));
        PetscCallAbort(
            communicator(),
            VecSet(ordered, 0.0));

        PetscInt sourceBegin = 0;
        PetscInt sourceEnd = 0;
        const PetscScalar *source = nullptr;
        PetscCallAbort(
            communicator(),
            VecGetOwnershipRange(
                currentOrderedVector,
                &sourceBegin,
                &sourceEnd));
        PetscCallAbort(
            communicator(),
            VecGetArrayRead(
                currentOrderedVector,
                &source));

        std::vector<PetscInt> destinationIndices(
            static_cast<std::size_t>(dofPerCell));
        std::vector<PetscScalar> values(
            static_cast<std::size_t>(dofPerCell));

        for (CellId cellId : ownedCells())
        {
            for (PetscInt component = 0; component < dofPerCell; ++component)
            {
                const PetscInt sourceGlobal =
                    globalDof(cellId, dofPerCell, component);

                if (sourceGlobal < sourceBegin || sourceGlobal >= sourceEnd)
                {
                    PetscCallAbort(
                        communicator(),
                        VecRestoreArrayRead(currentOrderedVector, &source));
                    PetscCallAbort(
                        communicator(),
                        VecDestroy(&ordered));
                    throw std::runtime_error(
                        "Owned StructuredGrid cell is outside Vec ownership range.");
                }

                destinationIndices[
                    static_cast<std::size_t>(component)] =
                    cellId * dofPerCell + component;
                values[static_cast<std::size_t>(component)] =
                    source[sourceGlobal - sourceBegin];
            }

            PetscCallAbort(
                communicator(),
                VecSetValues(
                    ordered,
                    dofPerCell,
                    destinationIndices.data(),
                    values.data(),
                    INSERT_VALUES));
        }

        PetscCallAbort(
            communicator(),
            VecRestoreArrayRead(currentOrderedVector, &source));
        PetscCallAbort(communicator(), VecAssemblyBegin(ordered));
        PetscCallAbort(communicator(), VecAssemblyEnd(ordered));
        return ordered;
    }

    /** @brief 结构网格的外部单元 id 即 Cartesian/current id。 */
    [[nodiscard]] PetscInt inputCellIndex(CellId cellId) const noexcept
    {
        return cellId;
    }

    [[nodiscard]] MPI_Comm communicator() const
    {
        return PetscObjectComm(
            reinterpret_cast<PetscObject>(
                grid_.dm(Indices::numPrimaryVariables)));
    }

    [[nodiscard]] PetscInt cellCount() const noexcept
    {
        const auto dim = grid_.dimensions();
        return static_cast<PetscInt>(dim[0]) *
               static_cast<PetscInt>(dim[1]) *
               static_cast<PetscInt>(dim[2]);
    }

    [[nodiscard]] bool isSetup() const noexcept
    {
        return grid_.isSetup();
    }

    [[nodiscard]] DofMap dofMap(PetscInt dof) const
    {
        return DofMap(*this, dof);
    }

    [[nodiscard]] std::vector<CellId> ownedCells() const
    {
        std::vector<CellId> result;
        result.reserve(grid_.cellIndices().size());
        for (int cell : grid_.cellIndices())
            result.push_back(static_cast<CellId>(cell));
        return result;
    }

    [[nodiscard]] std::vector<CellId> localSnapshotCells() const
    {
        const auto &info = layoutInfo_(Indices::numPrimaryVariables);
        std::vector<CellId> result;

        const std::size_t storageCount =
            static_cast<std::size_t>(info.ghost.xCount) *
            static_cast<std::size_t>(info.ghost.yCount) *
            static_cast<std::size_t>(info.ghost.zCount);
        result.reserve(storageCount);

        const auto insideOwned = [](int coordinate, int start, int count)
        {
            return coordinate >= start && coordinate < start + count;
        };

        for (int k = info.ghost.zStart;
             k < info.ghost.zStart + info.ghost.zCount;
             ++k)
        {
            for (int j = info.ghost.yStart;
                 j < info.ghost.yStart + info.ghost.yCount;
                 ++j)
            {
                for (int i = info.ghost.xStart;
                     i < info.ghost.xStart + info.ghost.xCount;
                     ++i)
                {
                    // PETSc 的 STAR local Vec 仍按完整长方体分配存储，但跨两个或
                    // 三个坐标方向的棱/角 ghost 不参与 GlobalToLocal 通信，值未定义。
                    // Natural 只缓存 owned 单元和实际通信的面 ghost；BOX stencil 则
                    // 可以使用完整 ghost 长方体。
                    if (info.stencilType == DMDA_STENCIL_STAR)
                    {
                        int outsideAxes = 0;
                        outsideAxes += !insideOwned(i, info.owned.xStart, info.owned.xCount);
                        outsideAxes += !insideOwned(j, info.owned.yStart, info.owned.yCount);
                        outsideAxes += !insideOwned(k, info.owned.zStart, info.owned.zCount);
                        if (outsideAxes > 1)
                            continue;
                    }

                    const int cell = grid_.ijkToGlobal(i, j, k);
                    result.push_back(static_cast<CellId>(cell));
                }
            }
        }

        return result;
    }

    [[nodiscard]] bool isLocal(CellId cellId) const noexcept
    {
        return grid_.onProcess(static_cast<int>(cellId));
    }

    [[nodiscard]] PetscInt localBlockIndex(
        CellId cellId,
        PetscInt dof) const
    {
        const auto &info = layoutInfo_(dof);
        const auto ijk = grid_.globalToIJK(
            static_cast<int>(cellId));

        const int i = ijk[0];
        const int j = ijk[1];
        const int k = ijk[2];

        if (i < info.ghost.xStart ||
            i >= info.ghost.xStart + info.ghost.xCount ||
            j < info.ghost.yStart ||
            j >= info.ghost.yStart + info.ghost.yCount ||
            k < info.ghost.zStart ||
            k >= info.ghost.zStart + info.ghost.zCount)
        {
            throw std::out_of_range(
                "Structured Natural cell is outside the DMDA ghost region.");
        }

        const PetscInt ii = i - info.ghost.xStart;
        const PetscInt jj = j - info.ghost.yStart;
        const PetscInt kk = k - info.ghost.zStart;

        return (kk * info.ghost.yCount + jj) *
                   info.ghost.xCount +
               ii;
    }

    [[nodiscard]] PetscInt localIndex(
        CellId cellId,
        PetscInt dof,
        PetscInt component) const
    {
        validateComponent_(dof, component);
        return localBlockIndex(cellId, dof) * dof + component;
    }

    [[nodiscard]] PetscInt globalDof(
        CellId cellId,
        PetscInt dof,
        PetscInt component) const
    {
        validateComponent_(dof, component);
        const auto &info = layoutInfo_(dof);
        const PetscInt local = localIndex(cellId, dof, component);
        PetscInt global = -1;

        PetscCallAbort(
            communicator(),
            ISLocalToGlobalMappingApply(
                info.localToGlobal,
                1,
                &local,
                &global));

        if (global < 0)
        {
            throw std::runtime_error(
                "Structured Natural DMDA local-to-global mapping returned an invalid DOF.");
        }

        return global;
    }

    [[nodiscard]] PetscInt ownedCellCount(PetscInt dof) const
    {
        return ownedDofCount(dof) / dof;
    }

    [[nodiscard]] PetscInt ghostCellCount(PetscInt dof) const
    {
        const auto &info = layoutInfo_(dof);
        const PetscInt localCells =
            static_cast<PetscInt>(info.ghost.xCount) *
            static_cast<PetscInt>(info.ghost.yCount) *
            static_cast<PetscInt>(info.ghost.zCount);
        return localCells - ownedCellCount(dof);
    }

    /**
     * @brief 返回 PETSc local Vec 的长方体 block 存储数量。
     *
     * 对 DMDA_STENCIL_STAR，该数量包含不会由 GlobalToLocal 通信的棱/角
     * ghost 存储槽；上层不得据此假设每个 block 都具有可读取的 ghost 值。
     */
    [[nodiscard]] PetscInt localCellCount(PetscInt dof) const
    {
        return ownedCellCount(dof) + ghostCellCount(dof);
    }

    [[nodiscard]] PetscInt ownedDofCount(PetscInt dof) const
    {
        const auto &info = layoutInfo_(dof);
        return info.ownedDofEnd - info.ownedDofBegin;
    }

    [[nodiscard]] PetscInt globalDofCount(PetscInt dof) const
    {
        return layoutInfo_(dof).globalDofCount;
    }

    [[nodiscard]] PetscInt ownedDofBegin(PetscInt dof) const
    {
        return layoutInfo_(dof).ownedDofBegin;
    }

    [[nodiscard]] PetscInt ownedBlockOrdinal(
        CellId cellId,
        PetscInt dof) const
    {
        const PetscInt global = globalDof(cellId, dof, 0);
        const PetscInt begin = ownedDofBegin(dof);
        if (global < begin || global >= begin + ownedDofCount(dof))
        {
            throw std::runtime_error(
                "Structured Natural requested owned block ordinal for a non-owned cell.");
        }
        return (global - begin) / dof;
    }

    [[nodiscard]] LocalReadView localReadView(
        PetscInt dof,
        Vec global) const
    {
        return LocalReadView(*this, dof, global);
    }

    [[nodiscard]] AssemblyView assemblyView() const
    {
        return AssemblyView(*this);
    }

    /**
     * @brief 返回 owned cell 的 DMDA 拓扑邻居，不访问岩石属性数组。
     */
    [[nodiscard]] std::vector<CellId> neighborCellIds(CellId cellId) const
    {
        const auto faces = grid_.faces(static_cast<int>(cellId));
        std::vector<CellId> result;
        result.reserve(faces.size());
        for (const auto &face : faces)
            result.push_back(static_cast<CellId>(face.neighbor()));
        return result;
    }

private:
    [[nodiscard]] const LayoutInfo &layoutInfo_(PetscInt dof) const
    {
        const auto found = layoutCache_.find(dof);
        if (found != layoutCache_.end())
            return *found->second;

        auto info = std::make_unique<LayoutInfo>();
        info->dof = dof;

        DM dm = grid_.dm(static_cast<int>(dof));

        PetscInt xs = 0, ys = 0, zs = 0;
        PetscInt xm = 0, ym = 0, zm = 0;
        PetscCallAbort(
            communicator(),
            DMDAGetCorners(dm, &xs, &ys, &zs, &xm, &ym, &zm));
        info->owned = {
            static_cast<int>(xs),
            static_cast<int>(ys),
            static_cast<int>(zs),
            static_cast<int>(xm),
            static_cast<int>(ym),
            static_cast<int>(zm)};

        PetscCallAbort(
            communicator(),
            DMDAGetGhostCorners(dm, &xs, &ys, &zs, &xm, &ym, &zm));
        info->ghost = {
            static_cast<int>(xs),
            static_cast<int>(ys),
            static_cast<int>(zs),
            static_cast<int>(xm),
            static_cast<int>(ym),
            static_cast<int>(zm)};

        PetscCallAbort(
            communicator(),
            DMDAGetInfo(
                dm, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                nullptr, nullptr, nullptr, nullptr, nullptr, &info->stencilType));

        PetscCallAbort(
            communicator(),
            DMGetLocalToGlobalMapping(dm, &info->localToGlobal));

        Vec probe = nullptr;
        PetscCallAbort(communicator(), DMCreateGlobalVector(dm, &probe));
        PetscCallAbort(
            communicator(),
            VecGetOwnershipRange(
                probe,
                &info->ownedDofBegin,
                &info->ownedDofEnd));
        PetscCallAbort(
            communicator(),
            VecGetSize(probe, &info->globalDofCount));
        PetscCallAbort(communicator(), VecDestroy(&probe));

        const auto result =
            layoutCache_.emplace(dof, std::move(info));
        return *result.first->second;
    }

    void validateAllCellsActive_() const
    {
        if (!grid_.isSetup())
        {
            throw std::logic_error(
                "StructuredGrid::setup() must be completed before Natural runtime construction.");
        }

        const auto owned = grid_.ownedRegion();
        const PetscInt expected =
            static_cast<PetscInt>(owned.xCount) *
            static_cast<PetscInt>(owned.yCount) *
            static_cast<PetscInt>(owned.zCount);

        if (static_cast<PetscInt>(grid_.cellIndices().size()) != expected)
        {
            throw std::invalid_argument(
                "Natural StructuredGrid runtime currently requires all DMDA cells to be active. "
                "Inactive-cell elimination must be handled by a reduced DOF layout, not by leaving zero rows in SNES.");
        }
    }

    static void validateComponent_(
        PetscInt dof,
        PetscInt component)
    {
        if (dof <= 0 || component < 0 || component >= dof)
        {
            throw std::out_of_range(
                "Natural StructuredGrid DOF component is out of range.");
        }
    }

    Grid &grid_;
    mutable std::unordered_map<
        PetscInt,
        std::unique_ptr<LayoutInfo>>
        layoutCache_;
};

} // namespace MPMC
