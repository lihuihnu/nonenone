/**
 * @file structuredgrid_legacy_geometry.hpp
 * @brief 旧 StructuredGrid<Tag> geometry Vec/raw-array API 的隔离兼容实现。
 */
#pragma once

#include <array>
#include <cstddef>

#include <petscdm.h>
#include <petscsys.h>
#include <petscvec.h>

#include <structuredgrid/structured_cartesian_geometry.hpp>

namespace MPMC::detail
{

/**
 * @brief 旧 `StructuredGrid<Tag>` façade 专用的 geometry Vec 兼容组件。
 *
 * 生产 `StructuredGridCore` 不继承该组件，因此不会创建 6/18-DOF geometry DMDA
 * 或 distance/areaNormal/center/volume Vec。该组件只为历史源码保留原字段与
 * `mapLocalArrays()` 后 raw-array 访问行为。
 */
class StructuredGridLegacyGeometryCompatibility
{
  public:
    Vec distance{nullptr};
    Vec areaNormal{nullptr};
    Vec center{nullptr};
    Vec volume{nullptr};

    std::array<double, 6> ***distanceArr{nullptr};
    std::array<double, 18> ***areaNormalArr{nullptr};
    std::array<double, 3> ***centerArr{nullptr};
    std::array<double, 1> ***volumeArr{nullptr};

  protected:
    template <class Grid>
    void initializeLegacyGeometry_(Grid &grid)
    {
        Vec distanceGlobal = grid.createGlobalVector(6);
        Vec areaNormalGlobal = grid.createGlobalVector(18);
        Vec centerGlobal = grid.createGlobalVector(3);
        Vec volumeGlobal = grid.createGlobalVector(1);

        auto distArray = grid.template vecGetArray<6>(distanceGlobal);
        auto areaArray = grid.template vecGetArray<18>(areaNormalGlobal);
        auto centerArray = grid.template vecGetArray<3>(centerGlobal);
        auto volumeArray = grid.template vecGetArray<1>(volumeGlobal);

        const auto centerX = structuredAxisCenters(grid.dx);
        const auto centerY = structuredAxisCenters(grid.dy);
        const auto centerZ = structuredAxisCenters(grid.dz);
        const auto owned = grid.ownedRegion();

        for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        {
            for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            {
                for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
                {
                    const std::array<double, 3> widths{
                        grid.dx[static_cast<std::size_t>(i)],
                        grid.dy[static_cast<std::size_t>(j)],
                        grid.dz[static_cast<std::size_t>(k)]};

                    distArray[k][j][i] = {
                        structuredFaceHalfDistance(widths, StructuredFaceOrder::XMinus),
                        structuredFaceHalfDistance(widths, StructuredFaceOrder::XPlus),
                        structuredFaceHalfDistance(widths, StructuredFaceOrder::YMinus),
                        structuredFaceHalfDistance(widths, StructuredFaceOrder::YPlus),
                        structuredFaceHalfDistance(widths, StructuredFaceOrder::ZMinus),
                        structuredFaceHalfDistance(widths, StructuredFaceOrder::ZPlus)};

                    centerArray[k][j][i] = {
                        centerX[static_cast<std::size_t>(i)],
                        centerY[static_cast<std::size_t>(j)],
                        centerZ[static_cast<std::size_t>(k)]};

                    areaArray[k][j][i] = {
                        structuredFaceNormalComponent(widths, StructuredFaceOrder::XMinus), 0.0, 0.0,
                        structuredFaceNormalComponent(widths, StructuredFaceOrder::XPlus), 0.0, 0.0,
                        0.0, structuredFaceNormalComponent(widths, StructuredFaceOrder::YMinus), 0.0,
                        0.0, structuredFaceNormalComponent(widths, StructuredFaceOrder::YPlus), 0.0,
                        0.0, 0.0, structuredFaceNormalComponent(widths, StructuredFaceOrder::ZMinus),
                        0.0, 0.0, structuredFaceNormalComponent(widths, StructuredFaceOrder::ZPlus)};

                    volumeArray[k][j][i][0] = structuredCellVolume(widths);
                }
            }
        }

        grid.template vecRestoreArray<6>(distanceGlobal, distArray);
        grid.template vecRestoreArray<18>(areaNormalGlobal, areaArray);
        grid.template vecRestoreArray<3>(centerGlobal, centerArray);
        grid.template vecRestoreArray<1>(volumeGlobal, volumeArray);

        distance = createLegacyLocal_(grid, 6, distanceGlobal);
        areaNormal = createLegacyLocal_(grid, 18, areaNormalGlobal);
        center = createLegacyLocal_(grid, 3, centerGlobal);
        volume = createLegacyLocal_(grid, 1, volumeGlobal);

        destroyLegacyVec_(distanceGlobal);
        destroyLegacyVec_(areaNormalGlobal);
        destroyLegacyVec_(centerGlobal);
        destroyLegacyVec_(volumeGlobal);
    }

    template <class Grid>
    void mapLegacyGeometryArrays_(Grid &grid)
    {
        distanceArr = grid.template vecGetArray<6>(distance);
        areaNormalArr = grid.template vecGetArray<18>(areaNormal);
        centerArr = grid.template vecGetArray<3>(center);
        volumeArr = grid.template vecGetArray<1>(volume);
    }

    template <class Grid>
    void restoreLegacyGeometryArrays_(Grid &grid) noexcept
    {
        if (distanceArr != nullptr)
            grid.template vecRestoreArray<6>(distance, distanceArr);
        if (areaNormalArr != nullptr)
            grid.template vecRestoreArray<18>(areaNormal, areaNormalArr);
        if (centerArr != nullptr)
            grid.template vecRestoreArray<3>(center, centerArr);
        if (volumeArr != nullptr)
            grid.template vecRestoreArray<1>(volume, volumeArr);
    }

    void destroyLegacyGeometry_() noexcept
    {
        destroyLegacyVec_(distance);
        destroyLegacyVec_(areaNormal);
        destroyLegacyVec_(center);
        destroyLegacyVec_(volume);
    }

  private:
    template <class Grid>
    [[nodiscard]] static Vec createLegacyLocal_(
        const Grid &grid,
        int dof,
        Vec global)
    {
        Vec local = nullptr;
        PetscCallAbort(PETSC_COMM_SELF,
                       DMCreateLocalVector(grid.dm(dof), &local));
        grid.globalToLocal(dof, global, local);
        return local;
    }

    static void destroyLegacyVec_(Vec &vec) noexcept
    {
        if (vec != nullptr)
            PetscCallAbort(PETSC_COMM_SELF, VecDestroy(&vec));
    }
};

} // namespace MPMC::detail
