/**
 * @file cpgrid_vtk_writer.hpp
 * @brief 将 CpGrid 单元结果写为 VTK 可视化文件。
 */
#pragma once

/*
 * Optional VTK backend.
 *
 * This header is intentionally NOT included by <output/output.hpp>.  Only
 * applications that actually request VTK output need VTK headers/libraries.
 */

#include <cpgrid/cpgrid.hpp>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <petscsys.h>

#include <cstddef>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace MPMC
{

enum class VtkCellOrdering
{
    CurrentId,
    InputIndex,
    StorageIndex = InputIndex
};

/**
 * @brief 面向稳定 CpGrid Mesh API 的 VTK 写出器。
 *
 * Geometry is inserted in original input-file order.  Field ordering is explicit
 * and is mapped to that geometry order, preventing ambiguity between canonical
 * input order and distributed current-id/DOF order.
 *
 * This class contains no MPI collectives.  For a full replicated CpGrid mesh,
 * call it on rank 0 after gathering the desired global field.
 */
template <class Grid>
class CpGridVtkWriter final
{
public:
    explicit CpGridVtkWriter(
        const Grid &grid)
        : grid_(grid),
          vtkGrid_(
              vtkSmartPointer<
                  vtkUnstructuredGrid>::New())
    {
        buildPoints_();
        buildCells_();
        addCellIdArrays_();
    }

    CpGridVtkWriter(
        const CpGridVtkWriter &) = delete;
    CpGridVtkWriter &operator=(
        const CpGridVtkWriter &) = delete;
    CpGridVtkWriter(
        CpGridVtkWriter &&) = delete;
    CpGridVtkWriter &operator=(
        CpGridVtkWriter &&) = delete;

    ~CpGridVtkWriter() = default;

    void addCellScalar(
        const std::string &name,
        const std::vector<double> &values,
        VtkCellOrdering ordering =
            VtkCellOrdering::InputIndex)
    {
        validateCellValueCount_(
            values.size());

        vtkNew<vtkDoubleArray> array;
        array->SetName(
            name.c_str());
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(
            static_cast<vtkIdType>(
                grid_.mesh().cellCount()));

        vtkIdType tuple = 0;

        for (const auto &cell :
             grid_.mesh().allCells())
        {
            const std::size_t source =
                fieldIndex_(
                    cell,
                    ordering);

            array->SetValue(
                tuple,
                values.at(source));

            ++tuple;
        }

        vtkGrid_
            ->GetCellData()
            ->AddArray(array);
    }

    /**
     * @brief 从 cell-major 扁平字段添加一个 VTK 分量。
     *
     * values is [cell0 c0..cN-1][cell1 c0..cN-1]... in the selected cell
     * ordering.
     */
    void addCellComponent(
        const std::string &name,
        const std::vector<double> &values,
        std::size_t componentsPerCell,
        std::size_t component,
        VtkCellOrdering ordering =
            VtkCellOrdering::InputIndex)
    {
        if (componentsPerCell == 0 ||
            component >= componentsPerCell)
        {
            throw std::invalid_argument(
                "Invalid VTK field component request.");
        }

        const std::size_t cells =
            grid_.mesh().cellCount();

        if (values.size() !=
            cells *
                componentsPerCell)
        {
            throw std::invalid_argument(
                "VTK cell-major field size does not match mesh cell count.");
        }

        vtkNew<vtkDoubleArray> array;
        array->SetName(
            name.c_str());
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(
            static_cast<vtkIdType>(
                cells));

        vtkIdType tuple = 0;

        for (const auto &cell :
             grid_.mesh().allCells())
        {
            const std::size_t sourceCell =
                fieldIndex_(
                    cell,
                    ordering);

            array->SetValue(
                tuple,
                values[
                    sourceCell *
                        componentsPerCell +
                    component]);

            ++tuple;
        }

        vtkGrid_
            ->GetCellData()
            ->AddArray(array);
    }

    void write(
        const std::string &filename) const
    {
        vtkNew<
            vtkXMLUnstructuredGridWriter>
            writer;

        writer->SetFileName(
            filename.c_str());

        writer->SetInputData(
            vtkGrid_);

        if (writer->Write() != 1)
        {
            throw std::runtime_error(
                "VTK XML unstructured-grid write failed.");
        }
    }

private:
    void buildPoints_()
    {
        vtkNew<vtkPoints> points;
        vtkNew<vtkIdTypeArray> nodeIds;

        nodeIds->SetName(
            "cpgrid_node_id");
        nodeIds->SetNumberOfComponents(1);

        nodeToVtk_.reserve(
            grid_.mesh().nodeCount());

        for (std::size_t index = 0;
             index < grid_.mesh().nodeCount();
             ++index)
        {
            const auto &node =
                grid_.mesh().node(index);

            const vtkIdType vtkId =
                points->InsertNextPoint(
                    node.x(),
                    node.y(),
                    node.z());

            const auto inserted =
                nodeToVtk_.emplace(
                    node.id(),
                    vtkId);

            if (!inserted.second)
            {
                throw std::runtime_error(
                    "CpGrid VTK output requires unique node ids.");
            }

            nodeIds->InsertNextValue(
                static_cast<vtkIdType>(
                    node.id()));
        }

        vtkGrid_->SetPoints(
            points);

        vtkGrid_
            ->GetPointData()
            ->AddArray(
                nodeIds);
    }

    void buildCells_()
    {
        for (const auto &cell :
             grid_.mesh().allCells())
        {
            std::set<vtkIdType>
                uniquePoints;

            vtkNew<vtkIdList>
                faceStream;

            for (const auto &face :
                 cell.faces())
            {
                faceStream->InsertNextId(
                    static_cast<vtkIdType>(
                        face.nodeCount()));

                for (const auto *node :
                     face.nodes())
                {
                    if (node == nullptr)
                    {
                        throw std::runtime_error(
                            "CpGrid VTK output encountered a null face node.");
                    }

                    const auto found =
                        nodeToVtk_.find(
                            node->id());

                    if (found ==
                        nodeToVtk_.end())
                    {
                        throw std::runtime_error(
                            "CpGrid VTK output cannot map a face node to VTK points.");
                    }

                    uniquePoints.insert(
                        found->second);

                    faceStream->InsertNextId(
                        found->second);
                }
            }

            std::vector<vtkIdType>
                pointIds(
                    uniquePoints.begin(),
                    uniquePoints.end());

            if (pointIds.empty() || cell.faceCount() == 0)
            {
                throw std::runtime_error(
                    "CpGrid VTK output encountered an empty polyhedron.");
            }

            vtkGrid_->InsertNextCell(
                VTK_POLYHEDRON,
                static_cast<vtkIdType>(
                    pointIds.size()),
                pointIds.data(),
                static_cast<vtkIdType>(
                    cell.faceCount()),
                faceStream->GetPointer(0));
        }
    }

    void addCellIdArrays_()
    {
        vtkNew<vtkIdTypeArray>
            currentIds;

        vtkNew<vtkIdTypeArray>
            inputIndices;

        currentIds->SetName(
            "cpgrid_current_id");

        inputIndices->SetName(
            "cpgrid_input_index");

        currentIds->SetNumberOfComponents(1);
        inputIndices->SetNumberOfComponents(1);

        for (const auto &cell :
             grid_.mesh().allCells())
        {
            currentIds->InsertNextValue(
                static_cast<vtkIdType>(
                    cell.id()));

            inputIndices->InsertNextValue(
                static_cast<vtkIdType>(
                    cell.inputIndex()));
        }

        vtkGrid_
            ->GetCellData()
            ->AddArray(
                currentIds);

        vtkGrid_
            ->GetCellData()
            ->AddArray(
                inputIndices);
    }

    void validateCellValueCount_(
        std::size_t size) const
    {
        if (size !=
            grid_.mesh().cellCount())
        {
            throw std::invalid_argument(
                "VTK scalar field size does not match mesh cell count.");
        }
    }

    [[nodiscard]] static std::size_t
    fieldIndex_(
        const typename Grid::Cell &cell,
        VtkCellOrdering ordering)
    {
        switch (ordering)
        {
        case VtkCellOrdering::CurrentId:
            return static_cast<std::size_t>(
                cell.id());

        case VtkCellOrdering::InputIndex:
            return static_cast<std::size_t>(
                cell.inputIndex());
        }

        throw std::logic_error(
            "Unknown VTK cell ordering.");
    }

    const Grid &grid_;
    vtkSmartPointer<
        vtkUnstructuredGrid>
        vtkGrid_;

    std::unordered_map<
        PetscInt,
        vtkIdType>
        nodeToVtk_;
};

} // namespace MPMC
