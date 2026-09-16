/**
 * @file aqueous_co2.hpp
 * @brief 水相 CO2 溶解度、Henry 关系及质量分数换算。
 */
#pragma once

#include <ad/Math.hpp>
#include <common/math.hpp>
#include <common/units.hpp>

#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 水相 CO2 的 Henry 定律热力学关系。
 *
 * 平衡条件写成 `f_CO2,w = x_CO2,w H(T,p_w,m_s)`，其中 `x_CO2,w`
 * 是水相 CO2 摩尔分数。Henry 常数采用压力修正
 * `H = H_ref(T) f_salt exp[Vbar_CO2 (p_w-p_sat)/(R T)]`。
 * `f_salt` 表示盐析修正，`Vbar_CO2` 为 CO2 在水中的偏摩尔体积。
 *
 * 该类只负责上述热力学关联，不知道网格、PETSc 或组分方程索引。
 */
class AqueousCO2Model final
{
public:
    AqueousCO2Model() = default;

    AqueousCO2Model(
        double temperature,
        double co2MolarMass,
        double waterMolarMass,
        double salinityMolality = 0.0)
        : temperature_(temperature),
          co2MolarMass_(co2MolarMass),
          waterMolarMass_(waterMolarMass),
          salinityMolality_(salinityMolality)
    {
        validate();
    }

    /** @brief 纯水在当前温度下的饱和蒸汽压 p_sat [Pa]。 */
    [[nodiscard]] double waterSaturationPressure() const
    {
        constexpr double Tc = 647.096;
        constexpr double Pc = 22.064e6;
        const double tau = 1.0 - temperature_ / Tc;
        const double a1 = -7.85951783;
        const double a2 = 1.84408259;
        const double a3 = -11.7866497;
        const double a4 = 22.6807411;
        const double a5 = -15.9618719;
        const double a6 = 1.80122502;
        const double exponent = (Tc / temperature_) *
            (a1 * tau + a2 * std::pow(tau, 1.5) + a3 * std::pow(tau, 3.0)
             + a4 * std::pow(tau, 3.5) + a5 * std::pow(tau, 4.0)
             + a6 * std::pow(tau, 7.5));
        return Pc * std::exp(exponent);
    }

    /** @brief CO2 在水中的偏摩尔体积 Vbar_CO2 [m^3/mol]。 */
    [[nodiscard]] double co2PartialMolarVolume() const
    {
        const double celsius = temperature_ - 273.15;
        const double cm3PerMol = 37.51 - 0.09585 * celsius
            + 8.740e-4 * celsius * celsius
            - 5.044e-7 * celsius * celsius * celsius;
        return cm3PerMol * 1.0e-6;
    }

    /** @brief 当前温度、纯水基准压力下的 Henry 常数 H_ref [Pa]。 */
    [[nodiscard]] double referenceHenryConstant() const
    {
        constexpr double Tc = 647.096;
        const double tr = temperature_ / Tc;
        if (!(tr > 0.0 && tr < 1.0))
            throw std::runtime_error("Aqueous CO2 Henry correlation requires 0 < T/Tc(water) < 1.");
        const double tau = 1.0 - tr;
        constexpr double A = -8.55445;
        constexpr double B = 4.01195;
        constexpr double C = 9.52345;
        const double exponent = A / tr
            + B * std::pow(tau, 0.355) / tr
            + C * std::pow(tr, -0.41) * std::exp(tau);
        return waterSaturationPressure() * std::exp(exponent);
    }

    /** @brief 盐度对 Henry 常数的无量纲盐析修正因子。 */
    [[nodiscard]] double saltingFactor() const
    {
        const double celsius = temperature_ - 273.15;
        const double ksalt =
            0.11572 - 6.0293e-4 * celsius
            + 3.5817e-6 * celsius * celsius
            - 3.7772e-9 * celsius * celsius * celsius;
        return std::pow(10.0, ksalt * salinityMolality_);
    }

    /** @brief 压力和盐度修正后的 Henry 常数 H(T,p_w,m_s) [Pa]。 */
    template <class Scalar>
    [[nodiscard]] Scalar henryConstant(const Scalar &waterPressure) const
    {
        constexpr double R = units::gasConstant;
        const double href = referenceHenryConstant();
        const double fsalt = saltingFactor();
        const double vbar = co2PartialMolarVolume();
        const double psat = waterSaturationPressure();
        return href * fsalt * exp(vbar * (waterPressure - psat) / (R * temperature_));
    }

    /**
     * @brief 将水相 CO2 摩尔分数转换为质量分数。
     *
     * `X_CO2 = x M_CO2 / [(1-x) M_H2O + x M_CO2]`。
     */
    template <class Scalar>
    [[nodiscard]] Scalar massFraction(const Scalar &co2MoleFraction) const
    {
        const Scalar numerator = co2MoleFraction * co2MolarMass_;
        const Scalar denominator =
            (1.0 - co2MoleFraction) * waterMolarMass_ + numerator;
        return numerator / denominator;
    }

    /** @brief 水相 CO2 逸度 `f_CO2,w = x_CO2,w H` [Pa]。 */
    template <class Scalar>
    [[nodiscard]] Scalar fugacity(
        const Scalar &waterPressure,
        const Scalar &co2MoleFraction) const
    {
        return co2MoleFraction * henryConstant(waterPressure);
    }

private:
    void validate() const
    {
        if (!(temperature_ > 0.0) || !std::isfinite(temperature_))
            throw std::invalid_argument("Aqueous CO2 temperature must be finite and positive.");
        if (!(co2MolarMass_ > 0.0) || !(waterMolarMass_ > 0.0))
            throw std::invalid_argument("Aqueous CO2 molar masses must be positive.");
        if (!std::isfinite(salinityMolality_) || salinityMolality_ < 0.0)
            throw std::invalid_argument("Salinity molality must be finite and nonnegative.");
    }

    double temperature_{298.15};
    double co2MolarMass_{0.04401};
    double waterMolarMass_{0.01801528};
    double salinityMolality_{0.0};
};

} // namespace MPMC
