/**
 * @file structuredgrid_views.hpp
 * @brief StructuredGrid 的受控 PETSc 岩石属性访问与装配期数组生命周期封装。
 */
#pragma once

#include <array>
#include <stdexcept>
#include <utility>

#include <petscsys.h>

namespace MPMC
{

/**
 * @brief StructuredGrid 全局岩石属性的可写 RAII 视图。
 *
 * 视图只负责当前 rank 拥有区域的全局 Vec 数组映射生命周期；调用方仍按
 * StructuredGrid::ownedRegion() 决定写入范围。析构时按照旧调用路径的顺序
 * 恢复渗透率和孔隙度数组映射。
 */
template <class Grid>
class StructuredGridRockWriteView final
{
  public:
    using VectorHandle = decltype(
        std::declval<const Grid &>().permeabilityVectorHandle());

    explicit StructuredGridRockWriteView(Grid &grid)
        : grid_(grid),
          permeabilityVector_(grid_.permeabilityVectorHandle()),
          porosityVector_(grid_.porosityVectorHandle())
    {
        if (!grid_.hasRockProperties())
            throw std::logic_error(
                "StructuredGrid rock properties must be initialized before write access.");

        permeability_ = grid_.template vecGetArray<3>(permeabilityVector_);
        porosity_ = grid_.template vecGetArray<1>(porosityVector_);
    }

    StructuredGridRockWriteView(const StructuredGridRockWriteView &) = delete;
    StructuredGridRockWriteView &operator=(const StructuredGridRockWriteView &) = delete;
    StructuredGridRockWriteView(StructuredGridRockWriteView &&) = delete;
    StructuredGridRockWriteView &operator=(StructuredGridRockWriteView &&) = delete;

    ~StructuredGridRockWriteView() noexcept
    {
        if (permeability_ != nullptr)
            grid_.template vecRestoreArray<3>(
                permeabilityVector_, permeability_);
        if (porosity_ != nullptr)
            grid_.template vecRestoreArray<1>(
                porosityVector_, porosity_);
    }

    [[nodiscard]] std::array<double, 3> &permeability(
        int i,
        int j,
        int k) noexcept
    {
        return permeability_[k][j][i];
    }

    [[nodiscard]] double &porosity(int i, int j, int k) noexcept
    {
        return porosity_[k][j][i][0];
    }

  private:
    Grid &grid_;
    VectorHandle permeabilityVector_{};
    VectorHandle porosityVector_{};
    std::array<double, 3> ***permeability_{nullptr};
    std::array<double, 1> ***porosity_{nullptr};
};

/**
 * @brief StructuredGrid 全局岩石属性的只读 RAII 视图。
 *
 * 该视图用于报告和诊断等只读取 owned global Vec 的场景，避免调用方直接
 * 管理 DMDAVecGetArrayRead()/DMDAVecRestoreArrayRead() 配对。
 */
template <class Grid>
class StructuredGridRockReadView final
{
  public:
    explicit StructuredGridRockReadView(const Grid &grid)
        : grid_(grid)
    {
        if (!grid_.hasRockProperties())
            throw std::logic_error(
                "StructuredGrid rock properties must be initialized before read access.");

        // 保持旧 grid report 的映射与恢复顺序：porosity -> permeability。
        porosity_ = grid_.template vecGetArrayRead<1>(
            grid_.porosityVectorHandle());
        permeability_ = grid_.template vecGetArrayRead<3>(
            grid_.permeabilityVectorHandle());
    }

    StructuredGridRockReadView(const StructuredGridRockReadView &) = delete;
    StructuredGridRockReadView &operator=(const StructuredGridRockReadView &) = delete;
    StructuredGridRockReadView(StructuredGridRockReadView &&) = delete;
    StructuredGridRockReadView &operator=(StructuredGridRockReadView &&) = delete;

    ~StructuredGridRockReadView() noexcept
    {
        if (porosity_ != nullptr)
            grid_.template vecRestoreArrayRead<1>(
                grid_.porosityVectorHandle(), porosity_);
        if (permeability_ != nullptr)
            grid_.template vecRestoreArrayRead<3>(
                grid_.permeabilityVectorHandle(), permeability_);
    }

    [[nodiscard]] std::array<double, 3> permeability(
        int i,
        int j,
        int k) const noexcept
    {
        return permeability_[k][j][i];
    }

    [[nodiscard]] double porosity(int i, int j, int k) const noexcept
    {
        return porosity_[k][j][i][0];
    }

  private:
    const Grid &grid_;
    const std::array<double, 3> ***permeability_{nullptr};
    const std::array<double, 1> ***porosity_{nullptr};
};

/**
 * @brief StructuredGrid 装配期 ghosted 岩石数组访问的 RAII 生命周期守卫。
 *
 * G7 后生产 Core 的 Cartesian 几何量按需计算；该类型只接管 ghosted
 * permeability/porosity 的 mapLocalArrays()/unmapLocalArrays() 严格配对。
 * 旧 `StructuredGrid<Tag>` façade 仍保留完整几何数组手工接口。
 */
template <class Grid>
class StructuredGridAssemblyAccess final
{
  public:
    explicit StructuredGridAssemblyAccess(Grid &grid)
        : grid_(grid)
    {
        grid_.mapLocalArrays();
        active_ = true;
    }

    StructuredGridAssemblyAccess(const StructuredGridAssemblyAccess &) = delete;
    StructuredGridAssemblyAccess &operator=(const StructuredGridAssemblyAccess &) = delete;
    StructuredGridAssemblyAccess(StructuredGridAssemblyAccess &&) = delete;
    StructuredGridAssemblyAccess &operator=(StructuredGridAssemblyAccess &&) = delete;

    ~StructuredGridAssemblyAccess() noexcept
    {
        if (!active_)
            return;

        try
        {
            grid_.unmapLocalArrays();
        }
        catch (...)
        {
            MPI_Abort(grid_.communicator(), PETSC_ERR_LIB);
        }
    }

  private:
    Grid &grid_;
    bool active_{false};
};

} // namespace MPMC
