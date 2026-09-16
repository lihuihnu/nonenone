/**
 * @file structuredgrid_components.hpp
 * @brief StructuredGrid 的内部几何、PETSc layout 与场存储组件。
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <petscsys.h>
#include <petscvec.h>

#include <structuredgrid/dmdacontainer.hpp>
#include <structuredgrid/structuredgridindex.hpp>
#include <structuredgrid/structured_cartesian_geometry.hpp>

namespace MPMC::detail
{

/**
 * @brief StructuredGrid 的纯几何定义与输入状态。
 *
 * 该组件不拥有 PETSc 资源，只保存全局尺寸、轴向单元宽度、自然编号器和
 * setup() 前的 active 输入，并集中负责纯几何输入校验。
 */
class StructuredGridGeometryState
{
  public:
    StructuredGridGeometryState(int nx, int ny, int nz)
        : nx_(nx), ny_(ny), nz_(nz), index_(nx, ny, nz)
    {
    }

    double Lx_{0.0};
    double Ly_{0.0};
    double Lz_{0.0};

    int nx_{0};
    int ny_{0};
    int nz_{0};

    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    std::vector<double> activeVec;

    [[nodiscard]] std::array<double, 3> physicalExtent() const noexcept
    {
        return {Lx_, Ly_, Lz_};
    }

    [[nodiscard]] StructuredGridIndex::Dimensions dimensions() const noexcept
    {
        return index_.dimensions();
    }

    [[nodiscard]] StructuredGridIndex::Coordinates
    globalToIJK(int globalIndex) const
    {
        return index_.globalToIJK(globalIndex);
    }

    [[nodiscard]] int ijkToGlobal(int i, int j, int k) const
    {
        return index_.ijkToGlobal(i, j, k);
    }

    [[nodiscard]] bool contains(int i, int j, int k) const noexcept
    {
        return index_.contains(i, j, k);
    }

    [[nodiscard]] std::array<double, 3>
    cellSize(const std::array<int, 3> &ijk) const
    {
        const auto [i, j, k] = ijk;
        return {dx.at(static_cast<std::size_t>(i)),
                dy.at(static_cast<std::size_t>(j)),
                dz.at(static_cast<std::size_t>(k))};
    }

  protected:
    StructuredGridIndex index_;

    void validateGridCounts_() const
    {
        if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0)
            throw std::invalid_argument(
                "StructuredGrid 三个方向的网格数量必须均大于 0");
    }

    static void validatePositiveLength_(double value, const char *name)
    {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument(std::string("StructuredGrid ") + name +
                                        " 必须为有限正数");
    }

    template <class AxisInput>
    static double initializeAxis_(const AxisInput &input,
                                  int count,
                                  std::vector<double> &output,
                                  const char *axisName)
    {
        using InputType = std::decay_t<AxisInput>;

        if constexpr (std::is_same_v<InputType, std::vector<double>>)
        {
            if (input.size() != static_cast<std::size_t>(count))
                throw std::invalid_argument(
                    std::string("StructuredGrid ") + axisName +
                    " 方向尺寸数组长度与网格数量不一致");

            for (double value : input)
                validatePositiveLength_(value, axisName);

            output = input;
            const double extent =
                std::accumulate(output.begin(), output.end(), 0.0);
            validatePositiveLength_(extent, axisName);
            return extent;
        }
        else if constexpr (std::is_arithmetic_v<InputType>)
        {
            const double value = static_cast<double>(input);
            validatePositiveLength_(value, axisName);
            output.assign(static_cast<std::size_t>(count), value);
            const double extent = value * count;
            validatePositiveLength_(extent, axisName);
            return extent;
        }
        else
        {
            static_assert(alwaysFalse_<InputType>,
                          "StructuredGrid 每个方向尺寸必须为算术标量或 std::vector<double>");
        }
    }

    void validateGeometry_()
    {
        validateGridCounts_();
        if (dx.size() != static_cast<std::size_t>(nx_) ||
            dy.size() != static_cast<std::size_t>(ny_) ||
            dz.size() != static_cast<std::size_t>(nz_))
            throw std::logic_error("StructuredGrid 单元尺寸数组尚未正确初始化");

        const auto validateAxis = [](const std::vector<double> &widths,
                                     const char *axisName) {
            for (double width : widths)
                validatePositiveLength_(width, axisName);
            const double extent =
                std::accumulate(widths.begin(), widths.end(), 0.0);
            validatePositiveLength_(extent, axisName);
            return extent;
        };

        Lx_ = validateAxis(dx, "x");
        Ly_ = validateAxis(dy, "y");
        Lz_ = validateAxis(dz, "z");
    }

    void initializeActivity_()
    {
        const auto totalCells = static_cast<std::size_t>(nx_) *
                                static_cast<std::size_t>(ny_) *
                                static_cast<std::size_t>(nz_);

        if (activeVec.empty())
        {
            activeVec.assign(totalCells, 1.0);
            return;
        }

        if (activeVec.size() != totalCells)
            throw std::invalid_argument(
                "StructuredGrid::activeVec 尺寸必须等于 nx * ny * nz");

        for (double &value : activeVec)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "StructuredGrid::activeVec 只能包含有限数值");
            value = value == 0.0 ? 0.0 : 1.0;
        }
    }

    void cacheVerticalCenters_()
    {
        verticalCenters_ = structuredAxisCenters(dz);
    }

    [[nodiscard]] double verticalCenter_(int k) const
    {
        return verticalCenters_.at(static_cast<std::size_t>(k));
    }

    std::vector<double> verticalCenters_;

  private:
    template <class>
    static constexpr bool alwaysFalse_{false};
};

/**
 * @brief StructuredGrid 的 DMDA 布局、并行区域与局部活动单元状态。
 */
class StructuredGridPetscLayoutState
{
  public:
    StructuredGridPetscLayoutState(int nx,
                                   int ny,
                                   int nz,
                                   MPI_Comm comm)
        : dmCont_(nx, ny, nz, 1, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, comm)
    {
    }

    int xl_{0};
    int yl_{0};
    int zl_{0};
    int nxl_{0};
    int nyl_{0};
    int nzl_{0};

    int gxl_{0};
    int gyl_{0};
    int gzl_{0};
    int gnxl_{0};
    int gnyl_{0};
    int gnzl_{0};

    void registerLayout(int dof)
    {
        dmCont_.registerLayout(dof);
    }

    [[nodiscard]] Vec createGlobalVector(int dof) const
    {
        return dmCont_.layout(dof).createGlobalVector();
    }

    [[nodiscard]] Vec borrowLocalVector(int dof) const
    {
        return dmCont_.layout(dof).borrowLocalVector();
    }

    [[nodiscard]] Vec borrowLocalVector(int dof, Vec global) const
    {
        return dmCont_.layout(dof).borrowLocalVector(global);
    }

    void restoreLocalVector(int dof, Vec &local) const
    {
        dmCont_.layout(dof).restoreLocalVector(local);
    }

    void globalToLocal(int dof, Vec global, Vec local) const
    {
        dmCont_.layout(dof).globalToLocal(global, local);
    }

    template <int dof>
    [[nodiscard]] std::array<double, dof> ***vecGetArray(Vec &vec) const
    {
        std::array<double, dof> ***arr = nullptr;
        PetscCallAbort(PETSC_COMM_SELF,
                       DMDAVecGetArray(dmCont_.layout(dof).dm(), vec, &arr));
        return arr;
    }

    template <std::size_t dof>
    void vecRestoreArray(Vec &vec, std::array<double, dof> ***&arr) const
    {
        PetscCallAbort(PETSC_COMM_SELF,
                       DMDAVecRestoreArray(
                           dmCont_.layout(static_cast<int>(dof)).dm(), vec, &arr));
        arr = nullptr;
    }

    template <int dof>
    [[nodiscard]] const std::array<double, dof> ***
    vecGetArrayRead(Vec vec) const
    {
        const std::array<double, dof> ***arr = nullptr;
        PetscCallAbort(PETSC_COMM_SELF,
                       DMDAVecGetArrayRead(dmCont_.layout(dof).dm(), vec, &arr));
        return arr;
    }

    template <std::size_t dof>
    void vecRestoreArrayRead(
        Vec vec, const std::array<double, dof> ***&arr) const
    {
        PetscCallAbort(
            PETSC_COMM_SELF,
            DMDAVecRestoreArrayRead(
                dmCont_.layout(static_cast<int>(dof)).dm(), vec, &arr));
        arr = nullptr;
    }

    [[nodiscard]] const std::vector<int> &cellIndices() const noexcept
    {
        return cellIndices_;
    }

    [[nodiscard]] DMDARegion ownedRegion() const noexcept
    {
        return DMDARegion{xl_, yl_, zl_, nxl_, nyl_, nzl_};
    }

    [[nodiscard]] bool isSetup() const noexcept
    {
        return setupComplete_;
    }

    [[nodiscard]] DM dm(int dof) const
    {
        return dmCont_.layout(dof).dm();
    }

    [[nodiscard]] Mat createMatrix(int dof) const
    {
        return dmCont_.layout(dof).createMatrix();
    }

    [[nodiscard]] MPI_Comm communicator() const noexcept
    {
        return dmCont_.communicator();
    }

  protected:
    [[nodiscard]] Vec createPersistentLocal_(int dof, Vec global) const
    {
        Vec local = nullptr;
        DM dmHandle = dmCont_.layout(dof).dm();
        PetscCallAbort(PETSC_COMM_SELF, DMCreateLocalVector(dmHandle, &local));
        dmCont_.layout(dof).globalToLocal(global, local);
        return local;
    }

    DMDAContainer dmCont_;
    std::vector<int> cellIndices_;
    std::vector<std::uint8_t> localActiveMask_;
    bool setupComplete_{false};
    bool arraysMapped_{false};
};

/**
 * @brief StructuredGrid 为 active ghost 同步保留的唯一几何相关 PETSc Vec。
 *
 * 正交 Cartesian 网格的半程距离、面面积、单元中心差和体积均可由轴宽
 * `dx/dy/dz` 按需计算，因此生产 Core 不再长期存储这些可推导几何场。
 */
struct StructuredGridActivityStorage
{
    Vec active{nullptr};
};

/**
 * @brief StructuredGrid 岩石属性全局 Vec 与装配期 ghosted 临时访问状态。
 */
struct StructuredGridRockStorage
{
    Vec porosityVector{nullptr};
    Vec permeabilityVector{nullptr};
    Vec localPermeabilityVector{nullptr};
    Vec localPorosityVector{nullptr};

    std::array<double, 3> ***permArr{nullptr};
    std::array<double, 1> ***refPoroArr{nullptr};
};

} // namespace MPMC::detail
