/**
 * @file structuredgridindex.hpp
 * @brief 结构网格自然编号、DMDA 编号与局部编号转换。
 */
#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 结构化网格的一维/三维单元索引转换器。
 *
 * 每个 StructuredGrid 独立持有一个 StructuredGridIndex，因此不同尺寸的
 * 结构化网格可以在同一进程中安全共存，不依赖任何全局或静态网格状态。
 *
 * 线性索引采用 x 方向最快、y 方向其次、z 方向最慢的布局：
 * `globalIndex = k * nx * ny + j * nx + i`。
 */
class StructuredGridIndex final
{
  public:
    using Index = int;
    using Coordinates = std::array<Index, 3>;
    using Dimensions = std::array<Index, 3>;

    /**
     * @brief 构造结构化网格索引器。
     *
     * @param[in] nx x 方向单元数量。
     * @param[in] ny y 方向单元数量。
     * @param[in] nz z 方向单元数量。
     *
     * @throws std::invalid_argument 当任一方向单元数量不为正数时抛出。
     * @throws std::overflow_error 当总单元数无法由 Index 表示时抛出。
     */
    StructuredGridIndex(Index nx, Index ny, Index nz)
        : nx_(nx), ny_(ny), nz_(nz)
    {
        validateDimensions_(nx_, ny_, nz_);
        planeSize_ = checkedProduct_(nx_, ny_);
        cellCount_ = checkedProduct_(planeSize_, nz_);
    }

    StructuredGridIndex() = delete;

    [[nodiscard]] constexpr Dimensions dimensions() const noexcept
    {
        return {nx_, ny_, nz_};
    }

    [[nodiscard]] constexpr Index nx() const noexcept
    {
        return nx_;
    }

    [[nodiscard]] constexpr Index ny() const noexcept
    {
        return ny_;
    }

    [[nodiscard]] constexpr Index nz() const noexcept
    {
        return nz_;
    }

    [[nodiscard]] constexpr Index planeSize() const noexcept
    {
        return planeSize_;
    }

    [[nodiscard]] constexpr Index cellCount() const noexcept
    {
        return cellCount_;
    }

    [[nodiscard]] constexpr bool contains(Index globalIndex) const noexcept
    {
        return globalIndex >= 0 && globalIndex < cellCount_;
    }

    [[nodiscard]] constexpr bool contains(Index i,
                                          Index j,
                                          Index k) const noexcept
    {
        return i >= 0 && i < nx_ && j >= 0 && j < ny_ && k >= 0 && k < nz_;
    }

    /**
     * @brief 将线性全局索引转换为 `(i,j,k)`。
     *
     * @throws std::out_of_range 当 globalIndex 越界时抛出。
     */
    [[nodiscard]] Coordinates globalToIJK(Index globalIndex) const
    {
        if (!contains(globalIndex))
        {
            throw std::out_of_range(
                "StructuredGridIndex::globalToIJK: global index is out of range.");
        }
        return globalToIJKUnchecked(globalIndex);
    }

    /**
     * @brief 将 `(i,j,k)` 转换为线性全局索引。
     *
     * @throws std::out_of_range 当任一方向索引越界时抛出。
     */
    [[nodiscard]] Index ijkToGlobal(Index i, Index j, Index k) const
    {
        if (!contains(i, j, k))
        {
            throw std::out_of_range(
                "StructuredGridIndex::ijkToGlobal: cell coordinates are out of range.");
        }
        return ijkToGlobalUnchecked(i, j, k);
    }

    /**
     * @brief 无边界检查地将线性索引转换为 `(i,j,k)`。
     *
     * @warning 仅在调用方已经证明索引有效时使用。
     */
    [[nodiscard]] constexpr Coordinates
    globalToIJKUnchecked(Index globalIndex) const noexcept
    {
        const Index k = globalIndex / planeSize_;
        const Index remainder = globalIndex - k * planeSize_;
        const Index j = remainder / nx_;
        const Index i = remainder - j * nx_;
        return {i, j, k};
    }

    /**
     * @brief 无边界检查地将 `(i,j,k)` 转换为线性索引。
     *
     * @warning 仅在调用方已经证明索引有效时使用。
     */
    [[nodiscard]] constexpr Index
    ijkToGlobalUnchecked(Index i, Index j, Index k) const noexcept
    {
        return k * planeSize_ + j * nx_ + i;
    }

  private:
    static void validateDimensions_(Index nx, Index ny, Index nz)
    {
        if (nx <= 0 || ny <= 0 || nz <= 0)
        {
            throw std::invalid_argument(
                "StructuredGridIndex requires positive grid dimensions.");
        }
    }

    [[nodiscard]] static Index checkedProduct_(Index lhs, Index rhs)
    {
        const auto product = static_cast<std::int64_t>(lhs) *
                             static_cast<std::int64_t>(rhs);
        if (product > static_cast<std::int64_t>(
                          std::numeric_limits<Index>::max()))
        {
            throw std::overflow_error(
                "StructuredGridIndex cell count exceeds the Index range.");
        }
        return static_cast<Index>(product);
    }

    Index nx_{0};
    Index ny_{0};
    Index nz_{0};
    Index planeSize_{0};
    Index cellCount_{0};
};

} // namespace MPMC
