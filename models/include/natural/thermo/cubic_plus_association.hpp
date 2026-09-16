/**
 * @file cubic_plus_association.hpp
 * @brief CPA 缔合 Helmholtz 项、位点分数和热力学修正。
 */
#pragma once

#include <array>
#include <cstddef>

namespace MPMC
{

/**
 * @brief CPA 缔合项使用的径向分布函数。
 *
 * `Simplified` 对应 sCPA：`g(eta)=1/(1-1.9 eta)`；
 * `CarnahanStarling` 对应硬球形式：`g(eta)=(1-eta/2)/(1-eta)^3`。
 */
enum class CpaRadialDistribution
{
    Simplified = 0,
    CarnahanStarling = 1
};

/**
 * @brief 常用 SRK-CPA 实现采用的标准 4C 水参数。
 *
 * 单位：`a0` [Pa·m6/mol2]、`b` [m3/mol]、`epsilon` [J/mol]。
 * 4C 水包含两个供体位点和两个受体位点。
 */
struct StandardCpaWater4C
{
    static constexpr double a0 = 0.12277;
    static constexpr double b = 1.4515e-5;
    static constexpr double c1 = 0.67359;
    static constexpr double epsilon = 16655.0;
    static constexpr double beta = 0.0692;
    static constexpr int donorSites = 2;
    static constexpr int acceptorSites = 2;
};

} // namespace MPMC
