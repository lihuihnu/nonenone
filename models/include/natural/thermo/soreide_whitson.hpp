/**
 * @file soreide_whitson.hpp
 * @brief Søreide–Whitson 水相修正参数与水相二元作用系数。
 */
#pragma once

#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief Søreide–Whitson PR 扩展使用的公开文献关联式。
 *
 * 参考：Soreide & Whitson, Fluid Phase Equilibria 77 (1992) 217–240；
 * Panfili et al., Computational Geosciences 28 (2024) 341–354。
 *
 * 物理：保留 Peng–Robinson 立方形式，只修改 H2O alpha，以及 Aqueous 相中的
 * H2O–溶质 BIP；盐度采用 NaCl molality [mol/kg H2O]。
 */
struct SoreideWhitsonCorrelations final
{
    /** @brief 原始 SW 水组分 alpha 关联式。 */
    template <class Scalar>
    [[nodiscard]] static Scalar waterAlpha(
        const Scalar &temperature,
        double waterCriticalTemperature,
        double salinityMolality)
    {
        if (!(waterCriticalTemperature > 0.0))
            throw std::invalid_argument("SW water critical temperature must be positive.");
        if (!(salinityMolality >= 0.0) || !std::isfinite(salinityMolality))
            throw std::invalid_argument("SW salinity must be finite and non-negative.");

        const Scalar tr = temperature / waterCriticalTemperature;
        const double saltFactor =
            1.0 - 0.0103 * std::pow(salinityMolality, 1.1);
        const Scalar sqrtAlpha =
            Scalar(1.0)
            + 0.4530 * (Scalar(1.0) - tr * saltFactor)
            + 0.0034 * (pow(tr, -3.0) - Scalar(1.0));
        return sqrtAlpha * sqrtAlpha;
    }



    /**
     * @brief 原始 SW 水相烃–H2O BIP 关联式。
     *
     * 约束：原拟合范围为 CH4 至 nC4；用于更重烃时属于外推。
     * 若重组分需要标定，应由算例提供专用 aqueous-BIP 回调。
     */
    [[nodiscard]] static double hydrocarbonAqueousBip(
        double temperature,
        double hydrocarbonCriticalTemperature,
        double acentricFactor,
        double salinityMolality)
    {
        if (!(temperature > 0.0) || !(hydrocarbonCriticalTemperature > 0.0))
            throw std::invalid_argument("SW hydrocarbon BIP temperatures must be positive.");
        if (!(acentricFactor > 0.0) || !std::isfinite(acentricFactor))
            throw std::invalid_argument("SW hydrocarbon acentric factor must be finite and positive.");
        if (!(salinityMolality >= 0.0) || !std::isfinite(salinityMolality))
            throw std::invalid_argument("SW salinity must be finite and non-negative.");

        const double tr = temperature / hydrocarbonCriticalTemperature;
        const double a0 = 1.1120 - 1.7369 * std::pow(acentricFactor, -0.1);
        const double a1 = 1.1001 + 0.8360 * acentricFactor;
        const double a2 = -0.15742 - 1.0988 * acentricFactor;

        // 文献：采用原始 SW 关联式中更正后的盐度系数。
        constexpr double alpha0 = 0.017407;
        constexpr double alpha1 = 0.033516;
        constexpr double alpha2 = 0.011478;

        return a0 * (1.0 + alpha0 * salinityMolality)
            + a1 * tr * (1.0 + alpha1 * salinityMolality)
            + a2 * tr * tr * (1.0 + alpha2 * salinityMolality);
    }

    /** @brief 原始 SW 水相 N2–H2O BIP 关联式。 */
    [[nodiscard]] static double nitrogenAqueousBip(
        double temperature,
        double nitrogenCriticalTemperature,
        double salinityMolality)
    {
        if (!(temperature > 0.0) || !(nitrogenCriticalTemperature > 0.0))
            throw std::invalid_argument("SW N2 BIP temperatures must be positive.");
        if (!(salinityMolality >= 0.0) || !std::isfinite(salinityMolality))
            throw std::invalid_argument("SW salinity must be finite and non-negative.");

        const double tr = temperature / nitrogenCriticalTemperature;
        const double salt = std::pow(salinityMolality, 0.75);
        return -1.70235 * (1.0 + 0.025587 * salt)
            + 0.44338 * (1.0 + 0.08126 * salt) * tr;
    }

    /**
     * @brief 原始 SW 水相 H2S–H2O BIP 关联式。
     *
     * 文献中的 H2S 表达式不含显式盐度项。
     */
    [[nodiscard]] static double hydrogenSulfideAqueousBip(
        double temperature,
        double hydrogenSulfideCriticalTemperature)
    {
        if (!(temperature > 0.0) || !(hydrogenSulfideCriticalTemperature > 0.0))
            throw std::invalid_argument("SW H2S BIP temperatures must be positive.");
        const double tr = temperature / hydrogenSulfideCriticalTemperature;
        return -0.20441 + 0.23426 * tr;
    }

    /**
     * @brief 原始 SW 水相 CO2–H2O BIP 关联式。
     *
     * `Tr` 为 CO2 约化温度；表达式对应 Panfili et al. (2024) Eq. (9)，
     * 其来源明确指向 Søreide–Whitson (1992) 原始 CO2 关联式。
     */
    [[nodiscard]] static double co2AqueousBip(
        double temperature,
        double co2CriticalTemperature,
        double salinityMolality)
    {
        if (!(temperature > 0.0) || !(co2CriticalTemperature > 0.0))
            throw std::invalid_argument("SW CO2 BIP temperatures must be positive.");
        if (!(salinityMolality >= 0.0) || !std::isfinite(salinityMolality))
            throw std::invalid_argument("SW salinity must be finite and non-negative.");

        const double tr = temperature / co2CriticalTemperature;
        return -0.31092 * (1.0 + 0.15587 * std::pow(salinityMolality, 0.7505))
            + 0.2358 * (1.0 + 0.17837 * std::pow(salinityMolality, 0.979)) * tr
            - 21.2566 * std::exp(-6.7222 * tr - salinityMolality);
    }

    /**
     * @brief Chabab et al. (2019) 改进的水相 CO2-H2O BIP 关联式。
     *
     * 原始 SW 参数在高 NaCl molality 下会系统性高估盐析效应。该关联式仅替换
     * aqueous CO2-H2O BIP，水 alpha、其他气体以及非水相 BIP 均保持原 SW
     * 定义。`Tr` 按论文固定为 `T / 304.13 K`，盐度为 mol/kg H2O。
     *
     * 参考：Chabab et al., International Journal of Greenhouse Gas Control 91
     * (2019) 102825, DOI 10.1016/j.ijggc.2019.102825。
     */
    [[nodiscard]] static double co2AqueousBipChabab2019(
        double temperature,
        double salinityMolality)
    {
        if (!(temperature > 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("Chabab SW CO2 BIP temperature must be finite and positive.");
        if (!(salinityMolality >= 0.0) || !std::isfinite(salinityMolality))
            throw std::invalid_argument("Chabab SW salinity must be finite and non-negative.");

        constexpr double co2CriticalTemperature = 304.13;
        constexpr double a = 0.43575155;
        constexpr double b = -0.05766906744;
        constexpr double c = 0.00826464849;
        constexpr double d = 0.00129539193;
        constexpr double e = -0.0016698848;
        constexpr double f = -0.47866096;
        const double tr = temperature / co2CriticalTemperature;
        return tr * (a + b * tr + c * tr * salinityMolality)
            + salinityMolality * salinityMolality * (d + e * tr)
            + f;
    }
};

} // namespace MPMC
