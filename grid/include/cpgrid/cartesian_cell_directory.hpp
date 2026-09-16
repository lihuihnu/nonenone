/**
 * @file cartesian_cell_directory.hpp
 * @brief 稀疏 Cartesian cell id 的紧凑只读目录。
 */
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace MPMC::detail
{

/**
 * @brief 保存 Cartesian id -> input index 映射，避免每 rank 持有全局 unordered_map。
 *
 * 常见的 GRDECL/MRST 数据中，input 顺序本身已经按 Cartesian id 递增，此时只
 * 保存一条 Cartesian-id 数组；若输入顺序不是单调的，则额外保存排序位置对应
 * 的 input index。查询使用二分，目录只读且内存连续。
 */
template <class Index>
class CompactCartesianCellDirectory final
{
    static_assert(std::is_integral_v<Index>, "Cartesian directory index must be integral.");

  public:
    void clear() noexcept
    {
        cartesianIds_.clear();
        inputIndexBySortedPosition_.clear();
    }

    void reset(
        const std::vector<Index> &cartesianByInput,
        const std::string &duplicateMessage)
    {
        clear();
        if (cartesianByInput.empty())
            return;

        bool inputOrderIsCartesianOrder = true;
        for (std::size_t input = 1; input < cartesianByInput.size(); ++input)
        {
            if (!(cartesianByInput[input - 1] < cartesianByInput[input]))
            {
                inputOrderIsCartesianOrder = false;
                break;
            }
        }

        if (inputOrderIsCartesianOrder)
        {
            cartesianIds_ = cartesianByInput;
            return;
        }

        std::vector<std::pair<Index, Index>> entries;
        entries.reserve(cartesianByInput.size());
        for (std::size_t input = 0; input < cartesianByInput.size(); ++input)
        {
            entries.emplace_back(
                cartesianByInput[input],
                static_cast<Index>(input));
        }
        std::sort(
            entries.begin(),
            entries.end(),
            [](const auto &left, const auto &right)
            {
                return left.first < right.first;
            });

        cartesianIds_.reserve(entries.size());
        inputIndexBySortedPosition_.reserve(entries.size());
        for (const auto &[cartesianId, inputIndex] : entries)
        {
            if (!cartesianIds_.empty() && cartesianIds_.back() == cartesianId)
                throw std::logic_error(duplicateMessage);
            cartesianIds_.push_back(cartesianId);
            inputIndexBySortedPosition_.push_back(inputIndex);
        }
    }

    [[nodiscard]] std::size_t size() const noexcept
    {
        return cartesianIds_.size();
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return cartesianIds_.empty();
    }

    /** @brief 常见的单调输入路径不需要额外 permutation 数组。 */
    [[nodiscard]] bool usesImplicitInputOrder() const noexcept
    {
        return inputIndexBySortedPosition_.empty();
    }

    [[nodiscard]] Index inputIndex(Index cartesianId) const
    {
        const auto found = std::lower_bound(
            cartesianIds_.begin(),
            cartesianIds_.end(),
            cartesianId);
        if (found == cartesianIds_.end() || *found != cartesianId)
            throw std::out_of_range("Cartesian cell id is not active in this Mesh.");

        const std::size_t sortedPosition =
            static_cast<std::size_t>(found - cartesianIds_.begin());
        if (inputIndexBySortedPosition_.empty())
            return static_cast<Index>(sortedPosition);
        return inputIndexBySortedPosition_[sortedPosition];
    }

    /**
     * @brief 仅为全局编号表/分发序列化临时重建 input-order Cartesian 数组。
     *
     * 生产求解不长期保存这条反向数组；调用方只在一次性诊断或 MPI 分发阶段
     * 承担 O(N) 临时内存。
     */
    [[nodiscard]] std::vector<Index> cartesianIdsInInputOrder() const
    {
        if (inputIndexBySortedPosition_.empty())
            return cartesianIds_;

        std::vector<Index> result(cartesianIds_.size());
        for (std::size_t sorted = 0; sorted < cartesianIds_.size(); ++sorted)
        {
            const Index input = inputIndexBySortedPosition_[sorted];
            if (input < 0 || static_cast<std::size_t>(input) >= result.size())
                throw std::logic_error("Cartesian directory contains an invalid input index.");
            result[static_cast<std::size_t>(input)] = cartesianIds_[sorted];
        }
        return result;
    }

  private:
    std::vector<Index> cartesianIds_;
    std::vector<Index> inputIndexBySortedPosition_;
};

} // namespace MPMC::detail
