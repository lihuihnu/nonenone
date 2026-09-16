/**
 * @file cpgrid_output_access.hpp
 * @brief CpGrid 单元输出数据的 owned/ghost 访问适配。
 */
#pragma once

#include <output/core/types.hpp>
#include <output/metrics/component_totals.hpp>
#include <output/petsc/owned_vector_view.hpp>

#include <cpgrid/cpgrid.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace MPMC
{

/**
 * @brief CpGrid PETSc 向量的只读输出适配器。
 *
 * 只读取 owned 条目：ghost 值不参与全局总量，采样值也只接受唯一 owner rank
 * 的贡献，从而避免 MPI 重复计数。
 */
template <class Grid>
class CpGridOutputAccess final
{
public:
    explicit CpGridOutputAccess(
        const Grid &grid)
        : grid_(grid)
    {
    }

    [[nodiscard]] std::vector<double>
    sample(
        Vec vector,
        PetscInt dofPerCell,
        const std::vector<OutputCellSelection> &cells,
        const std::vector<PetscInt> &components) const
    {
        validateVectorLayout_(
            vector,
            dofPerCell);

        if (components.empty())
        {
            return {};
        }

        for (PetscInt component : components)
        {
            if (component < 0 ||
                component >= dofPerCell)
            {
                throw std::out_of_range(
                    "Requested output component is outside the field layout.");
            }
        }

        const auto &map =
            grid_.dofMap(
                dofPerCell);

        OutputPetscOwnedReadView view(
            vector);

        const std::size_t valueCount =
            cells.size() *
            components.size();

        if (valueCount >
            static_cast<std::size_t>(
                std::numeric_limits<int>::max()))
        {
            throw std::overflow_error(
                "Output sample is too large for MPI_Allreduce count.");
        }

        std::vector<double>
            localValues(
                valueCount,
                0.0);

        std::vector<int>
            localOwners(
                cells.size(),
                0);

        for (std::size_t sampleIndex = 0;
             sampleIndex < cells.size();
             ++sampleIndex)
        {
            const PetscInt currentId =
                resolveCurrentCellId_(
                    cells[sampleIndex]);

            const PetscInt first =
                map.firstOwnedCellId();

            const PetscInt end =
                first +
                map.ownedCellCount();

            if (currentId < first ||
                currentId >= end)
            {
                continue;
            }

            localOwners[sampleIndex] = 1;

            const PetscInt localCell =
                currentId -
                first;

            for (std::size_t componentIndex = 0;
                 componentIndex < components.size();
                 ++componentIndex)
            {
                const PetscInt localDof =
                    localCell *
                        dofPerCell +
                    components[componentIndex];

                const std::size_t outputIndex =
                    sampleIndex *
                        components.size() +
                    componentIndex;

                localValues[outputIndex] =
                    PetscRealPart(
                        view.data()[
                            localDof]);
            }
        }

        std::vector<double>
            globalValues(
                valueCount,
                0.0);

        std::vector<int>
            globalOwners(
                cells.size(),
                0);

        const MPI_Comm comm =
            grid_.mesh().communicator();

        if (valueCount > 0)
        {
            PetscCallMPIAbort(
                comm,
                MPI_Allreduce(
                    localValues.data(),
                    globalValues.data(),
                    static_cast<int>(
                        valueCount),
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm));
        }

        if (!cells.empty())
        {
            PetscCallMPIAbort(
                comm,
                MPI_Allreduce(
                    localOwners.data(),
                    globalOwners.data(),
                    static_cast<int>(
                        cells.size()),
                    MPI_INT,
                    MPI_SUM,
                    comm));
        }

        for (std::size_t index = 0;
             index < globalOwners.size();
             ++index)
        {
            if (globalOwners[index] != 1)
            {
                throw std::runtime_error(
                    "Output cell sample does not have exactly one MPI owner.");
            }
        }

        return globalValues;
    }

    template <std::size_t NumComponents>
    [[nodiscard]] ComponentTotals<NumComponents>
    componentTotals(
        Vec vector) const
    {
        static_assert(
            NumComponents > 0,
            "Output component total requires at least one component.");

        constexpr PetscInt dofPerCell =
            static_cast<PetscInt>(
                NumComponents);

        validateVectorLayout_(
            vector,
            dofPerCell);

        const auto &map =
            grid_.dofMap(
                dofPerCell);

        OutputPetscOwnedReadView view(
            vector);

        ComponentTotals<NumComponents>
            local;

        for (PetscInt cell = 0;
             cell < map.ownedCellCount();
             ++cell)
        {
            for (std::size_t component = 0;
                 component < NumComponents;
                 ++component)
            {
                const PetscInt localDof =
                    cell *
                        dofPerCell +
                    static_cast<PetscInt>(
                        component);

                local.component[component] +=
                    PetscRealPart(
                        view.data()[
                            localDof]);
            }
        }

        ComponentTotals<NumComponents>
            global;

        const MPI_Comm comm =
            grid_.mesh().communicator();

        PetscCallMPIAbort(
            comm,
            MPI_Allreduce(
                local.component.data(),
                global.component.data(),
                static_cast<int>(
                    NumComponents),
                MPI_DOUBLE,
                MPI_SUM,
                comm));

        return global;
    }

private:
    void validateVectorLayout_(
        Vec vector,
        PetscInt dofPerCell) const
    {
        if (vector == nullptr)
        {
            throw std::invalid_argument(
                "CpGrid output requires a valid PETSc Vec.");
        }

        if (dofPerCell <= 0)
        {
            throw std::invalid_argument(
                "Output DOF per cell must be positive.");
        }

        if (!grid_.hasLayout(
                dofPerCell))
        {
            throw std::invalid_argument(
                "Requested CpGrid output layout has not been registered.");
        }

        PetscInt globalSize = 0;
        PetscInt localSize = 0;

        PetscCallAbort(
            grid_.mesh().communicator(),
            VecGetSize(
                vector,
                &globalSize));

        PetscCallAbort(
            grid_.mesh().communicator(),
            VecGetLocalSize(
                vector,
                &localSize));

        const auto &map =
            grid_.dofMap(
                dofPerCell);

        if (globalSize !=
                map.globalDofCount() ||
            localSize !=
                map.ownedDofCount())
        {
            throw std::invalid_argument(
                "PETSc Vec layout does not match the requested CpGrid DOF layout.");
        }
    }

    [[nodiscard]] PetscInt
    resolveCurrentCellId_(
        const OutputCellSelection &selection) const
    {
        if (selection.id < 0)
        {
            throw std::out_of_range(
                "Output cell id cannot be negative.");
        }

        switch (selection.space)
        {
        case OutputCellIdSpace::CurrentId:
        {
            const PetscInt currentId = static_cast<PetscInt>(selection.id);
            (void)grid_.mesh().inputIndexFromCurrentId(currentId);
            return currentId;
        }

        case OutputCellIdSpace::InputIndex:
            return grid_.mesh().currentIdFromInputIndex(
                static_cast<std::size_t>(selection.id));

        case OutputCellIdSpace::CartesianId:
            return grid_.mesh().currentIdFromCartesianId(
                static_cast<PetscInt>(selection.id));
        }

        throw std::logic_error(
            "Unknown output cell identifier space.");
    }

    const Grid &grid_;
};

} // namespace MPMC
