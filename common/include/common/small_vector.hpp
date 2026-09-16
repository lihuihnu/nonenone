/**
 * @file small_vector.hpp
 * @brief 面向小维度数据的轻量连续容器。
 */
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 带固定内联缓冲区的小型向量。
 *
 * 当元素数量不超过 `InlineCapacity` 时，数据直接保存在对象内部；
 * 超过容量时退化为 `std::vector`。
 *
 * 本实现不缓存裸 `dataPtr_`，从而避免 move 后
 * 指针失效、自移动赋值和动态/内联存储切换时的悬空指针风险。
 *
 * @tparam ValueType 元素类型。
 * @tparam InlineCapacity 内联缓冲区容量。
 */
template <class ValueType, std::size_t InlineCapacity>
class SmallVector final
{
  public:
    SmallVector() = default;

    explicit SmallVector(std::size_t size)
    {
        resize_(size);
    }

    SmallVector(
        std::size_t size,
        const ValueType &value)
    {
        resize_(size);
        std::fill(
            data(),
            data() + size_,
            value);
    }

    SmallVector(const SmallVector &) = default;
    SmallVector &operator=(const SmallVector &) = default;

    SmallVector(SmallVector &&other) noexcept(
        std::is_nothrow_move_constructible_v<ValueType>)
    {
        moveFrom_(std::move(other));
    }

    SmallVector &operator=(SmallVector &&other) noexcept(
        std::is_nothrow_move_assignable_v<ValueType>)
    {
        if (this != &other)
        {
            moveFrom_(std::move(other));
        }

        return *this;
    }

    ~SmallVector() = default;

    [[nodiscard]] std::size_t size() const noexcept
    {
        return size_;
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return size_ == 0;
    }

    [[nodiscard]] ValueType *data() noexcept
    {
        return usesInlineStorage_()
                   ? inlineBuffer_.data()
                   : heapBuffer_.data();
    }

    [[nodiscard]] const ValueType *data() const noexcept
    {
        return usesInlineStorage_()
                   ? inlineBuffer_.data()
                   : heapBuffer_.data();
    }

    [[nodiscard]] ValueType &operator[](
        std::size_t index) noexcept
    {
        return data()[index];
    }

    [[nodiscard]] const ValueType &operator[](
        std::size_t index) const noexcept
    {
        return data()[index];
    }

    [[nodiscard]] ValueType &at(
        std::size_t index)
    {
        if (index >= size_)
        {
            throw std::out_of_range(
                "SmallVector index is out of range.");
        }

        return (*this)[index];
    }

    [[nodiscard]] const ValueType &at(
        std::size_t index) const
    {
        if (index >= size_)
        {
            throw std::out_of_range(
                "SmallVector index is out of range.");
        }

        return (*this)[index];
    }

    [[nodiscard]] ValueType *begin() noexcept
    {
        return data();
    }

    [[nodiscard]] ValueType *end() noexcept
    {
        return size_ == 0
                   ? data()
                   : data() + size_;
    }

    [[nodiscard]] const ValueType *begin() const noexcept
    {
        return data();
    }

    [[nodiscard]] const ValueType *end() const noexcept
    {
        return size_ == 0
                   ? data()
                   : data() + size_;
    }

  private:
    [[nodiscard]] bool usesInlineStorage_() const noexcept
    {
        return size_ <= InlineCapacity;
    }

    void resize_(std::size_t size)
    {
        size_ = size;

        if (!usesInlineStorage_())
        {
            heapBuffer_.resize(size_);
        }
    }

    void moveFrom_(SmallVector &&other)
    {
        size_ = other.size_;

        if (other.usesInlineStorage_())
        {
            inlineBuffer_ =
                std::move(
                    other.inlineBuffer_);

            heapBuffer_.clear();
        }
        else
        {
            heapBuffer_ =
                std::move(
                    other.heapBuffer_);
        }

        other.size_ = 0;
        other.heapBuffer_.clear();
    }

    std::array<ValueType, InlineCapacity>
        inlineBuffer_{};

    std::vector<ValueType>
        heapBuffer_;

    std::size_t size_{0};
};

} // namespace MPMC
