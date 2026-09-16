/**
 * @file filter_range.hpp
 * @brief 对容器视图进行轻量条件过滤的范围适配器。
 */
#pragma once

#include <iterator>
#include <utility>

namespace MPMC
{

/**
 * @brief C++17 轻量过滤迭代器。
 *
 * 用于在不复制 Polyhedron 的前提下遍历当前 MPI rank 拥有的单元。
 */
template <class Iterator, class Predicate>
class FilterIterator
{
  public:
    using iterator_category =
        std::forward_iterator_tag;
    using value_type =
        typename std::iterator_traits<Iterator>::value_type;
    using difference_type =
        typename std::iterator_traits<Iterator>::difference_type;
    using pointer =
        typename std::iterator_traits<Iterator>::pointer;
    using reference =
        typename std::iterator_traits<Iterator>::reference;

    FilterIterator(
        Iterator current,
        Iterator end,
        Predicate predicate)
        : current_(current),
          end_(end),
          predicate_(std::move(predicate))
    {
        advance_();
    }

    reference operator*() const
    {
        return *current_;
    }

    pointer operator->() const
    {
        return std::addressof(*current_);
    }

    FilterIterator &operator++()
    {
        ++current_;
        advance_();
        return *this;
    }

    FilterIterator operator++(int)
    {
        FilterIterator copy(*this);
        ++(*this);
        return copy;
    }

    friend bool operator==(
        const FilterIterator &lhs,
        const FilterIterator &rhs)
    {
        return lhs.current_ == rhs.current_;
    }

    friend bool operator!=(
        const FilterIterator &lhs,
        const FilterIterator &rhs)
    {
        return !(lhs == rhs);
    }

  private:
    void advance_()
    {
        while (current_ != end_ &&
               !predicate_(*current_))
        {
            ++current_;
        }
    }

    Iterator current_;
    Iterator end_;
    Predicate predicate_;
};

/**
 * @brief 由 begin/end 与谓词组成的可 range-for 过滤视图。
 */
template <class Iterator, class Predicate>
class FilterRange
{
  public:
    FilterRange(
        Iterator begin,
        Iterator end,
        Predicate predicate)
        : begin_(begin, end, predicate),
          end_(end, end, predicate)
    {
    }

    [[nodiscard]] auto begin() const
    {
        return begin_;
    }

    [[nodiscard]] auto end() const
    {
        return end_;
    }

  private:
    FilterIterator<Iterator, Predicate> begin_;
    FilterIterator<Iterator, Predicate> end_;
};

template <class Iterator, class Predicate>
FilterRange(
    Iterator,
    Iterator,
    Predicate)
    -> FilterRange<Iterator, Predicate>;

} // namespace MPMC
