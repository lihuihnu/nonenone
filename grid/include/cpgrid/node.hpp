/**
 * @file node.hpp
 * @brief CpGrid 网格节点数据结构。
 */
#pragma once

#include <cpgrid/point.hpp>

#include <petscsys.h>

namespace MPMC
{

/**
 * @brief 带稳定输入编号的网格节点。
 */
class Node final : public Point
{
  public:
    Node() = default;

    Node(
        double x,
        double y,
        double z,
        PetscInt id) noexcept
        : Point(x, y, z),
          id_(id)
    {
    }

    [[nodiscard]] PetscInt id() const noexcept
    {
        return id_;
    }

  private:
    PetscInt id_{0};
};

} // namespace MPMC
