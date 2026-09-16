/**
 * @file dmdacontainer.hpp
 * @brief StructuredGrid 的 PETSc DMDA 对象所有权与访问封装。
 */
#pragma once

#include <stdexcept>
#include <cstddef>
#include <map>

#include <structuredgrid/dmdaelem.hpp>

namespace MPMC
{

/**
 * @brief 缓存同一结构化网格拓扑上不同 dof 的 DMDA 布局。
 *
 * 第一个注册的 dof 创建 canonical DMDA；后续 dof 通过
 * DMDACreateCompatibleDMDA() 从 canonical 布局派生，从而保证所有物理场、
 * 几何场和辅助场具有完全一致的 MPI 分区与 ghost 拓扑。
 */
class DMDAContainer final
{
  public:
    using Dof = int;
    using Container = std::map<Dof, DMDAElem>;
    using Size = Container::size_type;

    DMDAContainer(int nx,
                  int ny,
                  int nz,
                  int stencilWidth = 1,
                  DMBoundaryType boundaryType = DM_BOUNDARY_NONE,
                  DMDAStencilType stencilType = DMDA_STENCIL_STAR,
                  MPI_Comm comm = PETSC_COMM_WORLD);

    DMDAContainer(const DMDAContainer &) = delete;
    DMDAContainer &operator=(const DMDAContainer &) = delete;
    DMDAContainer(DMDAContainer &&) noexcept = default;
    DMDAContainer &operator=(DMDAContainer &&) noexcept = default;

    /**
     * @brief 注册指定 dof 的数据布局。
     *
     * 已注册 dof 时保持幂等并直接返回。
     *
     * @note 首次注册以及创建兼容 DMDA 都是 MPI collective；所有 rank 必须以
     *       相同顺序执行首次注册操作。
     */
    void registerLayout(Dof dof);

    [[nodiscard]] bool contains(Dof dof) const noexcept
    {
        return layouts_.find(dof) != layouts_.end();
    }

    [[nodiscard]] Size size() const noexcept
    {
        return layouts_.size();
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return layouts_.empty();
    }

    /**
     * @brief 获取指定 dof 的 DMDA 布局。
     *
     * @throws std::out_of_range 当 dof 尚未注册时抛出。
     */
    [[nodiscard]] DMDAElem &layout(Dof dof)
    {
        return layouts_.at(dof);
    }

    [[nodiscard]] const DMDAElem &layout(Dof dof) const
    {
        return layouts_.at(dof);
    }

    [[nodiscard]] DMDADimensions requestedDimensions() const noexcept
    {
        return {nx_, ny_, nz_};
    }

    /** @brief 返回所有已注册 DMDA 共同使用的 MPI communicator。 */
    [[nodiscard]] MPI_Comm communicator() const noexcept
    {
        return comm_;
    }

  private:
    void validatePrimaryLayout_(const DMDAElem &layout) const;

    int nx_{0};
    int ny_{0};
    int nz_{0};
    int stencilWidth_{1};
    DMBoundaryType boundaryType_{DM_BOUNDARY_NONE};
    DMDAStencilType stencilType_{DMDA_STENCIL_STAR};
    MPI_Comm comm_{PETSC_COMM_WORLD};
    Dof canonicalDof_{0};
    Container layouts_;
};

} // namespace MPMC
