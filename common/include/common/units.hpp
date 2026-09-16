/**
 * @file units.hpp
 * @brief 油藏工程与 SI 单位之间的统一换算常量。
 */
#pragma once

namespace MPMC::units
{

/**
 * @brief MPMC 内部推荐使用 SI；以下常量仅用于经验关联式和输入换算。
 */
inline constexpr double metre = 1.0;
inline constexpr double second = 1.0;
inline constexpr double kilogram = 1.0;
inline constexpr double pascal = 1.0;

/** CODATA molar gas constant, J/(mol K). */
inline constexpr double gasConstant = 8.31446261815324;

inline constexpr double bar =
    100000.0 * pascal;

inline constexpr double poise =
    0.1 * pascal * second;

inline constexpr double centipoise =
    1.0e-2 * poise;

inline constexpr double pound =
    0.45359237 * kilogram;

inline constexpr double standardGravity =
    9.80665 * metre /
    (second * second);

inline constexpr double poundForce =
    pound * standardGravity;

inline constexpr double inch =
    0.0254 * metre;

inline constexpr double psi =
    poundForce /
    (inch * inch);

/**
 * @brief Rankine 数值转换为 Kelvin 数值时的比例系数。
 *
 * `T_K = T_R * kelvinPerRankine`。
 */
inline constexpr double kelvinPerRankine =
    5.0 / 9.0;

} // namespace MPMC::units
