/**
 * @file types.hpp
 * @brief CpGrid 公共类型、索引和轻量别名。
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string_view>

namespace MPMC
{

/** @brief 井的物理角色；储层方程中注入为正源、采出为负源。 */
enum class WellType
{
    Injector,
    Producer
};

/**
 * @brief 当前活动的井控制方程类型。
 *
 * `Bhp` constrains bottom-hole pressure. ReservoirTotalRate constrains total
 * in-situ volume; the other rate controls use surface-condition volume.
 */
enum class WellControl
{
    Bhp,
    TotalRate,
    ReservoirTotalRate,
    OilRate,
    GasRate,
    WaterRate
};

/** @brief 当控制残差属于流量方程时返回 true。 */
[[nodiscard]] constexpr bool isRateControl(WellControl control) noexcept
{
    return control != WellControl::Bhp;
}

/** @brief 返回文本/CSV 诊断使用的稳定大写标签。 */
[[nodiscard]] constexpr std::string_view wellTypeName(WellType type) noexcept
{
    switch (type)
    {
    case WellType::Injector: return "INJECTOR";
    case WellType::Producer: return "PRODUCER";
    }
    return "UNKNOWN";
}

/** @brief 返回文本/CSV 诊断使用的稳定控制标签。 */
[[nodiscard]] constexpr std::string_view wellControlName(WellControl control) noexcept
{
    switch (control)
    {
    case WellControl::Bhp:       return "BHP";
    case WellControl::TotalRate: return "RATE";
    case WellControl::ReservoirTotalRate: return "RES_RATE";
    case WellControl::OilRate:   return "OIL_RATE";
    case WellControl::GasRate:   return "GAS_RATE";
    case WellControl::WaterRate: return "WATER_RATE";
    }
    return "UNKNOWN";
}

/**
 * @brief 将非负工程流量目标转换为模拟器符号约定。
 *
 * Injection into the reservoir is positive; production from the reservoir is negative.
 */
[[nodiscard]] inline double signedRateTarget(WellType type, double targetMagnitude)
{
    if (!std::isfinite(targetMagnitude))
        throw std::invalid_argument("Well rate target must be finite.");
    const double magnitude = std::abs(targetMagnitude);
    return type == WellType::Producer ? -magnitude : magnitude;
}

/** @brief 将带符号模拟器流量转换为非负工程量。 */
[[nodiscard]] inline double rateMagnitude(WellType type, double signedRate) noexcept
{
    return type == WellType::Injector
        ? std::max(0.0, signedRate)
        : std::max(0.0, -signedRate);
}

} // namespace MPMC
