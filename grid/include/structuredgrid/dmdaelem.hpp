/**
 * @file dmdaelem.hpp
 * @brief StructuredGrid 单元、面和邻接关系辅助。
 */
#pragma once

#include <petscdmda.h>
#include <petscsys.h>

namespace MPMC
{

/**
 * @brief 三维 DMDA 的局部区域描述。
 */
struct DMDARegion
{
    int xStart{0};
    int yStart{0};
    int zStart{0};
    int xCount{0};
    int yCount{0};
    int zCount{0};
};

/**
 * @brief 三维 DMDA 的全局尺寸。
 */
struct DMDADimensions
{
    int nx{0};
    int ny{0};
    int nz{0};

    [[nodiscard]] friend constexpr bool
    operator==(const DMDADimensions &lhs, const DMDADimensions &rhs) noexcept
    {
        return lhs.nx == rhs.nx && lhs.ny == rhs.ny && lhs.nz == rhs.nz;
    }

    [[nodiscard]] friend constexpr bool
    operator!=(const DMDADimensions &lhs, const DMDADimensions &rhs) noexcept
    {
        return !(lhs == rhs);
    }
};

/**
 * @brief PETSc 三维 DMDA 的唯一资源拥有包装器。
 *
 * DMDAElem 独占一个 PETSc DM 句柄，不允许复制，只允许移动。这样底层 DM 的
 * 所有权在 C++ 层始终唯一，避免引入不必要的 PETSc 引用计数共享语义。
 *
 * @note 短期局部向量通过 borrowLocalVector() 从 PETSc 向量池借用，必须使用
 *       restoreLocalVector() 归还；长期局部向量应由调用层使用
 *       DMCreateLocalVector() 创建并通过 VecDestroy() 销毁。
 */
class DMDAElem final
{
  public:
    /**
     * @brief 创建并初始化一个三维 PETSc DMDA。
     *
     * @param[in] nx x 方向全局网格数量。
     * @param[in] ny y 方向全局网格数量。
     * @param[in] nz z 方向全局网格数量。
     * @param[in] dof 每个网格单元上的自由度数量。
     * @param[in] stencilWidth stencil/ghost 宽度。
     * @param[in] boundaryType 三个方向统一使用的边界类型。
     * @param[in] stencilType DMDA stencil 类型。
     * @param[in] comm MPI 通信器。
     */
    DMDAElem(int nx,
             int ny,
             int nz,
             int dof,
             int stencilWidth = 1,
             DMBoundaryType boundaryType = DM_BOUNDARY_NONE,
             DMDAStencilType stencilType = DMDA_STENCIL_STAR,
             MPI_Comm comm = PETSC_COMM_WORLD);

    /**
     * @brief 基于已有 DMDA 的并行布局创建不同 dof 的兼容 DMDA。
     *
     * @param[in] layoutSource 提供网格拓扑、并行划分和 stencil 的源 DMDA。
     * @param[in] dof 新 DMDA 的每单元自由度数量。
     *
     * @note 该操作内部使用 DMDACreateCompatibleDMDA()，因此属于 MPI collective。
     */
    DMDAElem(const DMDAElem &layoutSource, int dof);

    DMDAElem(const DMDAElem &) = delete;
    DMDAElem &operator=(const DMDAElem &) = delete;

    DMDAElem(DMDAElem &&other) noexcept;
    DMDAElem &operator=(DMDAElem &&other) noexcept;

    ~DMDAElem() noexcept;

    /**
     * @brief 获取每个单元的自由度数量。
     */
    [[nodiscard]] int dof() const noexcept
    {
        return dof_;
    }

    /**
     * @brief 获取非拥有型 PETSc DM 句柄。
     *
     * @warning 调用方不得对返回句柄调用 DMDestroy()。
     */
    [[nodiscard]] DM dm() const noexcept
    {
        return dm_;
    }

    /**
     * @brief 获取 DMDA 全局尺寸。
     */
    [[nodiscard]] DMDADimensions dimensions() const noexcept
    {
        return dimensions_;
    }

    /**
     * @brief 获取当前 MPI rank 拥有的非 ghost 区域。
     */
    [[nodiscard]] DMDARegion ownedRegion() const noexcept
    {
        return ownedRegion_;
    }

    /**
     * @brief 获取当前 MPI rank 含 ghost 的局部区域。
     */
    [[nodiscard]] DMDARegion ghostRegion() const noexcept
    {
        return ghostRegion_;
    }

    /**
     * @brief 创建拥有型全局 Vec。
     *
     * @return 新创建 Vec；调用方负责 VecDestroy()。
     */
    [[nodiscard]] Vec createGlobalVector() const;

    /**
     * @brief 从 DM 局部向量池借用一个短期 ghost Vec。
     *
     * @return 借用 Vec；必须通过 restoreLocalVector() 归还。
     */
    [[nodiscard]] Vec borrowLocalVector() const;

    /**
     * @brief 借用局部 Vec 并立即完成 global-to-local 同步。
     *
     * @param[in] global 输入全局 Vec。
     * @return 已同步 ghost 值的借用局部 Vec。
     */
    [[nodiscard]] Vec borrowLocalVector(Vec global) const;

    /**
     * @brief 将短期借用局部 Vec 归还给 PETSc DM。
     */
    void restoreLocalVector(Vec &local) const;

    /**
     * @brief 将全局 Vec 同步到已有局部 Vec。
     */
    void globalToLocal(Vec global, Vec local) const;

    /**
     * @brief 创建与当前 DMDA stencil 匹配的 PETSc Mat。
     *
     * @return 新创建矩阵；调用方负责 MatDestroy()。
     */
    [[nodiscard]] Mat createMatrix() const;

  private:
    void refreshMetadata_();
    void release_() noexcept;
    void reset_() noexcept;

    DM dm_{nullptr};
    int dof_{0};
    DMDADimensions dimensions_{};
    DMDARegion ownedRegion_{};
    DMDARegion ghostRegion_{};
};

} // namespace MPMC
