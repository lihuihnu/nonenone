/**
 * @file gridenums.hpp
 * @brief 结构网格方向、边界和索引枚举。
 */
#pragma once

#include <cstddef>
#include <cstdint>

namespace MPMC
{

/**
 * @brief 网格拓扑类型。
 *
 * 该枚举只描述网格的拓扑组织方式，不包含任何物理模型信息。
 */
enum class GridType : std::uint8_t
{
  Structured,     ///< 规则或逻辑规则的结构化网格。
  SemiStructured, ///< 半结构化网格。
  Unstructured    ///< 非结构化网格。
};

/**
 * @brief 结构化网格单元六个面的固定存储顺序。
 *
 * 为避免 `left/right/front/back/up/down` 与具体坐标系约定产生歧义，
 * 统一使用笛卡尔坐标轴及正负方向描述六个面。
 *
 * 该顺序与 StructuredGrid 中几何数组的布局保持一致：
 *
 * - 0: x 负方向面；
 * - 1: x 正方向面；
 * - 2: y 负方向面；
 * - 3: y 正方向面；
 * - 4: z 负方向面；
 * - 5: z 正方向面。
 *
 * `Count` 仅用于表示面数量，不代表实际网格面。
 */
enum class StructuredFaceOrder : std::uint8_t
{
  XMinus = 0,
  XPlus = 1,
  YMinus = 2,
  YPlus = 3,
  ZMinus = 4,
  ZPlus = 5,
  Count = 6
};

/**
 * @brief 判断给定枚举值是否表示有效的结构化网格面。
 *
 * @param[in] face 待检查的面枚举。
 * @return 若 face 对应六个实际面之一则返回 true，否则返回 false。
 */
[[nodiscard]] constexpr bool
isValidStructuredFace(StructuredFaceOrder face) noexcept
{
  return static_cast<std::uint8_t>(face) <
         static_cast<std::uint8_t>(StructuredFaceOrder::Count);
}

/**
 * @brief 将结构化网格面转换为几何数组索引。
 *
 * @param[in] face 结构化网格面。
 * @return 对应的零基索引；`Count` 对应 6。
 *
 * @note 调用方在索引六面数组前应确保 `isValidStructuredFace(face)` 为 true。
 */
[[nodiscard]] constexpr std::size_t
structuredFaceIndex(StructuredFaceOrder face) noexcept
{
  return static_cast<std::size_t>(face);
}

/**
 * @brief 获取与给定面相对的结构化网格面。
 *
 * @param[in] face 结构化网格面。
 * @return 相对面；若传入 `Count` 或无效枚举值，则返回 `Count`。
 */
[[nodiscard]] constexpr StructuredFaceOrder
oppositeStructuredFace(StructuredFaceOrder face) noexcept
{
  switch (face)
  {
  case StructuredFaceOrder::XMinus:
    return StructuredFaceOrder::XPlus;
  case StructuredFaceOrder::XPlus:
    return StructuredFaceOrder::XMinus;
  case StructuredFaceOrder::YMinus:
    return StructuredFaceOrder::YPlus;
  case StructuredFaceOrder::YPlus:
    return StructuredFaceOrder::YMinus;
  case StructuredFaceOrder::ZMinus:
    return StructuredFaceOrder::ZPlus;
  case StructuredFaceOrder::ZPlus:
    return StructuredFaceOrder::ZMinus;
  case StructuredFaceOrder::Count:
    return StructuredFaceOrder::Count;
  }

  return StructuredFaceOrder::Count;
}

/**
 * @brief 获取结构化网格面所属的笛卡尔坐标轴。
 *
 * @param[in] face 结构化网格面。
 * @return x、y、z 方向分别返回 0、1、2；无效面返回 -1。
 */
[[nodiscard]] constexpr int
structuredFaceAxis(StructuredFaceOrder face) noexcept
{
  switch (face)
  {
  case StructuredFaceOrder::XMinus:
  case StructuredFaceOrder::XPlus:
    return 0;
  case StructuredFaceOrder::YMinus:
  case StructuredFaceOrder::YPlus:
    return 1;
  case StructuredFaceOrder::ZMinus:
  case StructuredFaceOrder::ZPlus:
    return 2;
  case StructuredFaceOrder::Count:
    return -1;
  }

  return -1;
}

/**
 * @brief 获取结构化网格面相对于所属坐标轴的方向符号。
 *
 * @param[in] face 结构化网格面。
 * @return 负方向面返回 -1，正方向面返回 +1，无效面返回 0。
 */
[[nodiscard]] constexpr int
structuredFaceSign(StructuredFaceOrder face) noexcept
{
  switch (face)
  {
  case StructuredFaceOrder::XMinus:
  case StructuredFaceOrder::YMinus:
  case StructuredFaceOrder::ZMinus:
    return -1;
  case StructuredFaceOrder::XPlus:
  case StructuredFaceOrder::YPlus:
  case StructuredFaceOrder::ZPlus:
    return 1;
  case StructuredFaceOrder::Count:
    return 0;
  }

  return 0;
}

} // namespace MPMC