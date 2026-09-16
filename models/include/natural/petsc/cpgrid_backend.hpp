/**
 * @file cpgrid_backend.hpp
 * @brief Natural 运行时面向 CpGrid 的网格后端。
 */
#pragma once

#include <stdexcept>
#include <cpgrid/cpgrid.hpp>
#include <natural/petsc/grid_backend_common.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <cstddef>
#include <vector>

namespace MPMC
{

/**
 * @brief CpGrid 到 Natural PETSc runtime 的后端适配器。
 *
 * 后端只负责网格/PETSc 数据访问；不计算 EOS、相平衡、蓄积、井源项。
 */
template <class Indices, class Grid>
class NaturalCpGridBackend final
{
public:
    using CellId = PetscInt;
    using GridType = Grid;
    using DofMap = NaturalDofMapView<NaturalCpGridBackend>;

    class LocalReadView final
    {
    public:
        LocalReadView(
            const NaturalCpGridBackend &backend,
            PetscInt dofPerCell,
            Vec global)
            : backend_(backend),
              dofPerCell_(dofPerCell),
              local_(backend_.grid_.borrowLocalVector(
                  dofPerCell_,
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
                backend_.grid_.restoreLocalVector(dofPerCell_, local_);
        }

        [[nodiscard]] double value(
            CellId cellId,
            PetscInt component) const
        {
            const PetscInt localIndex =
                backend_.grid_.dofMap(dofPerCell_)
                    .localIndex(cellId, component);
            return PetscRealPart(array_[localIndex]);
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
        const NaturalCpGridBackend &backend_;
        PetscInt dofPerCell_{0};
        Vec local_{nullptr};
        const PetscScalar *array_{nullptr};
    };

    class AssemblyView final
    {
    public:
        explicit AssemblyView(
            const NaturalCpGridBackend &backend)
            : backend_(backend),
              rock_(backend_.grid_.rockLocalView())
        {
        }

        AssemblyView(const AssemblyView &) = delete;
        AssemblyView &operator=(const AssemblyView &) = delete;
        AssemblyView(AssemblyView &&) = delete;
        AssemblyView &operator=(AssemblyView &&) = delete;

        [[nodiscard]] double porosity(CellId cellId) const
        {
            return rock_.porosity(
                backend_.grid_.mesh().cellByCurrentId(cellId));
        }

        [[nodiscard]] double cellVolume(CellId cellId) const
        {
            return backend_.grid_.cellVolume(
                backend_.grid_.mesh().cellByCurrentId(cellId));
        }

        [[nodiscard]] std::vector<NaturalCellConnection>
        connections(CellId cellId) const
        {
            const auto &cell =
                backend_.grid_.mesh().cellByCurrentId(cellId);

            std::vector<NaturalCellConnection> result;
            result.reserve(cell.faces().size());

            for (const auto &face : cell.faces())
            {
                const auto *neighbor = face.neighborCell();
                if (neighbor == nullptr)
                    continue;

                result.push_back({
                    neighbor->id(),
                    backend_.grid_.transmissibility(face),
                    backend_.grid_.gravityTerm(face)});
            }

            return result;
        }

    private:
        const NaturalCpGridBackend &backend_;
        typename Grid::RockLocalView rock_;
    };

    explicit NaturalCpGridBackend(Grid &grid)
        : grid_(grid)
    {
    }

    [[nodiscard]] Grid &grid() noexcept { return grid_; }
    [[nodiscard]] const Grid &grid() const noexcept { return grid_; }

    void registerLayout(PetscInt dof)
    {
        grid_.registerLayout(dof);
    }

    /** @brief CpGrid 历史路径没有额外 AD 布局；显式保持 no-op。 */
    void registerAdAuxiliaryLayout(PetscInt) noexcept
    {
    }

    [[nodiscard]] Vec createGlobalVector(PetscInt dof) const
    {
        return grid_.createGlobalVector(dof);
    }

    /**
     * @brief 将求解器 Vec 复制为原始网格输入顺序用于外部输出。
     *
     * The solver keeps partition/current-id ordering internally.  This API is the
     * only ordering conversion exposed to the Natural runtime.
     */
    [[nodiscard]] Vec createInputOrderedCopy(
        PetscInt dof,
        Vec currentOrderedVector) const
    {
        return grid_.createInputOrderedCopy(
            dof,
            currentOrderedVector);
    }

    /** @brief 返回内部 current cell id 对应的原始零基网格文件行号。 */
    [[nodiscard]] PetscInt inputCellIndex(CellId currentCellId) const
    {
        return grid_.mesh().inputIndexFromCurrentId(currentCellId);
    }

    [[nodiscard]] MPI_Comm communicator() const noexcept
    {
        return grid_.mesh().communicator();
    }

    [[nodiscard]] PetscInt cellCount() const noexcept
    {
        return static_cast<PetscInt>(grid_.mesh().cellCount());
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
        result.reserve(
            static_cast<std::size_t>(
                grid_.dofMap(Indices::numPrimaryVariables)
                    .ownedCellCount()));

        const auto &owned = grid_.mesh().ownedCellIds();
        result.assign(owned.begin(), owned.end());

        return result;
    }

    [[nodiscard]] std::vector<CellId> localSnapshotCells() const
    {
        const auto &map =
            grid_.dofMap(
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables));

        std::vector<CellId> result(
            static_cast<std::size_t>(
                map.ownedCellCount() +
                map.ghostCellCount()),
            CellId(-1));

        const auto store = [&](const auto &cell)
        {
            result[static_cast<std::size_t>(
                map.localBlockIndex(cell))] =
                cell.id();
        };

        for (PetscInt ownedCellId : grid_.mesh().ownedCellIds())
            store(grid_.mesh().cellByCurrentId(ownedCellId));

        for (PetscInt ghostCellId : map.ghostCellIds())
            store(grid_.mesh().cellByCurrentId(ghostCellId));

        for (CellId cellId : result)
        {
            if (cellId < 0)
                throw std::runtime_error(
                    "CpGrid local snapshot cell table is incomplete.");
        }

        return result;
    }

    [[nodiscard]] bool isLocal(CellId cellId) const
    {
        return grid_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables))
            .isOwnedCellId(cellId);
    }

    [[nodiscard]] PetscInt localBlockIndex(
        CellId cellId,
        PetscInt dof) const
    {
        return grid_.dofMap(dof).localBlockIndex(cellId);
    }

    [[nodiscard]] PetscInt localIndex(
        CellId cellId,
        PetscInt dof,
        PetscInt component) const
    {
        return grid_.dofMap(dof).localIndex(cellId, component);
    }

    [[nodiscard]] PetscInt globalDof(
        CellId cellId,
        PetscInt dof,
        PetscInt component) const
    {
        return grid_.dofMap(dof).globalIndex(cellId, component);
    }

    [[nodiscard]] PetscInt ownedCellCount(PetscInt dof) const
    {
        return grid_.dofMap(dof).ownedCellCount();
    }

    [[nodiscard]] PetscInt ghostCellCount(PetscInt dof) const
    {
        return grid_.dofMap(dof).ghostCellCount();
    }

    [[nodiscard]] PetscInt localCellCount(PetscInt dof) const
    {
        return ownedCellCount(dof) + ghostCellCount(dof);
    }

    [[nodiscard]] PetscInt ownedDofCount(PetscInt dof) const
    {
        return grid_.dofMap(dof).ownedDofCount();
    }

    [[nodiscard]] PetscInt globalDofCount(PetscInt dof) const
    {
        return grid_.dofMap(dof).globalDofCount();
    }

    [[nodiscard]] PetscInt ownedDofBegin(PetscInt dof) const
    {
        return grid_.dofMap(dof).firstOwnedGlobalDof();
    }

    [[nodiscard]] PetscInt ownedBlockOrdinal(
        CellId cellId,
        PetscInt dof) const
    {
        const PetscInt first = ownedDofBegin(dof);
        const PetscInt global = globalDof(cellId, dof, 0);
        return (global - first) / dof;
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
     * @brief 返回 owned cell 的拓扑邻居，不访问岩石属性。
     *
     * Jacobian 预分配只需要连接关系，因此不应为了拓扑查询映射 K/poro。
     */
    [[nodiscard]] std::vector<CellId> neighborCellIds(CellId cellId) const
    {
        const auto &cell = grid_.mesh().cellByCurrentId(cellId);
        std::vector<CellId> result;
        result.reserve(cell.faces().size());
        for (const auto &face : cell.faces())
        {
            const auto *neighbor = face.neighborCell();
            if (neighbor != nullptr)
                result.push_back(neighbor->id());
        }
        return result;
    }

private:
    Grid &grid_;
};


} // namespace MPMC
