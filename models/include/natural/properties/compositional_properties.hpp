/**
 * @file compositional_properties.hpp
 * @brief 由 EOS 相组成计算密度、质量分数和黏度的组分物性。
 */
#pragma once

#include <common/units.hpp>
#include <natural/compositional_mixture.hpp>

#include <array>
#include <tuple>

namespace MPMC
{

/**
 * @brief EOS相密度与 Lohrenz-Bray-Clark (LBC) 风格黏度关联式。
 *
 * 组分数固定在编译期，组成量使用 `std::array`。SI 单位贯穿接口：
 * pressure [Pa]、temperature [K]、mass density [kg/m^3]、viscosity [Pa s]。
 */
template <class Indices>
class CompositionalPropertyModel
{
public:
    static constexpr int numComponents = Indices::numComponents;

    explicit CompositionalPropertyModel(const CompositionalMixture<Indices> &mixture)
        : mixture_(mixture)
    {
        initializeLbcConstants_();
    }

    CompositionalPropertyModel() = default;

    /** @brief EOS 摩尔密度 `rho_m = p/(R T Z)` [mol/m^3]. */
    template <class Scalar>
    [[nodiscard]] Scalar molarDensity(
        const Scalar &pressure,
        const Scalar &compressibility,
        const Scalar &temperature) const
    {
        return pressure / (units::gasConstant * temperature * compressibility);
    }

    /**
     * @brief 相质量密度 `rho = rho_m * M_mix` [kg/m^3].
     *
     * `M_mix = sum_i x_i M_i` 为当前相平均摩尔质量。
     */
    template <class Scalar>
    [[nodiscard]] Scalar densityFromMolarDensity(
        const std::array<Scalar, numComponents> &composition,
        const Scalar &molarDensityValue) const
    {
        Scalar mixtureMolarMass = 0.0;
        for (int i = 0; i < numComponents; ++i)
            mixtureMolarMass += composition[static_cast<std::size_t>(i)] * mixture_.molecularWeight(i);
        return molarDensityValue * mixtureMolarMass;
    }

    template <class Scalar>
    [[nodiscard]] Scalar density(
        const Scalar &pressure,
        const std::array<Scalar, numComponents> &composition,
        const Scalar &compressibility,
        const Scalar &temperature) const
    {
        return densityFromMolarDensity(
            composition, molarDensity(pressure, compressibility, temperature));
    }

    /**
     * @brief Kay 型混合伪临界参数 `(P_pc,T_pc,V_pc,M_mix)`.
     *
     * 每个量采用摩尔分数线性混合，例如 `T_pc=sum_i x_i T_ci`。
     */
    template <class Scalar>
    [[nodiscard]] std::tuple<Scalar, Scalar, Scalar, Scalar>
    pseudoCriticalProperties(const std::array<Scalar, numComponents> &composition) const
    {
        Scalar pseudoCriticalPressure = 0.0;
        Scalar pseudoCriticalTemperature = 0.0;
        Scalar pseudoCriticalVolume = 0.0;
        Scalar mixtureMolarMass = 0.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const auto x = composition[static_cast<std::size_t>(i)];
            pseudoCriticalPressure += mixture_.criticalPressure(i) * x;
            pseudoCriticalTemperature += mixture_.criticalTemperature(i) * x;
            pseudoCriticalVolume += mixture_.criticalVolume(i) * x;
            mixtureMolarMass += mixture_.molecularWeight(i) * x;
        }
        return {pseudoCriticalPressure, pseudoCriticalTemperature,
                pseudoCriticalVolume, mixtureMolarMass};
    }

    /**
     * @brief LBC 风格相黏度 [Pa s].
     *
     * 先按低压组分黏度求混合值，再用约化摩尔密度
     * `rho_r = V_pc rho_m` 计算高压密度修正。该函数只依赖 `(p,T,Z,x)`，
     * 不需要知道相标签；液/气差异已经包含在各自的 `Z` 与组成中。
     */
    template <class Scalar>
    [[nodiscard]] Scalar viscosityFromMolarDensity(
        const std::array<Scalar, numComponents> &composition,
        const Scalar &rhoMolar,
        const Scalar &temperature) const
    {
        constexpr double molarMassScale = 1000.0;

        Scalar numerator = 0.0;
        Scalar denominator = 0.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const double mwSqrt = lbcMolarMassSqrt_[idx];
            const Scalar reducedTemperature = temperature / mixture_.criticalTemperature(i);
            const double epsilon = lbcEpsilon_[idx];

            const Scalar componentViscosity = scalarValue(reducedTemperature) > 1.5
                ? 17.78e-5 * pow(4.58 * reducedTemperature - 1.67, 0.625) / epsilon
                : 34.0e-5 * pow(reducedTemperature, 0.94) / epsilon;
            numerator += composition[static_cast<std::size_t>(i)] * componentViscosity * mwSqrt;
            denominator += composition[static_cast<std::size_t>(i)] * mwSqrt;
        }
        const Scalar atmosphericViscosity = numerator / denominator;

        auto [pcMix, tcMix, vcMix, mwMix] = pseudoCriticalProperties(composition);
        const Scalar epsilonMix =
            5.4402 * pow(tcMix / units::kelvinPerRankine, 1.0 / 6.0) /
            (sqrt(molarMassScale * mwMix) * pow(pcMix / units::psi, 2.0 / 3.0) * units::centipoise);

        const Scalar reducedDensity = vcMix * rhoMolar;
        constexpr std::array<double, 6> c{
            0.1023, 0.023364, 0.058533, -0.040758, 0.0093324, -1.0e-4};
        // Integer powers are hot in cell-property evaluation.  Explicit
        // multiplication is both clearer and cheaper for primitive and AD scalars
        // than routing these exact powers through the generic pow() overload.
        const Scalar reducedDensity2 = reducedDensity * reducedDensity;
        const Scalar reducedDensity3 = reducedDensity2 * reducedDensity;
        const Scalar reducedDensity4 = reducedDensity2 * reducedDensity2;
        const Scalar polynomial =
            c[0] + c[1] * reducedDensity + c[2] * reducedDensity2
            + c[3] * reducedDensity3 + c[4] * reducedDensity4;
        const Scalar polynomial2 = polynomial * polynomial;
        const Scalar polynomial4 = polynomial2 * polynomial2;
        return atmosphericViscosity + (polynomial4 + c[5]) / epsilonMix;
    }

    template <class Scalar>
    [[nodiscard]] Scalar viscosity(
        const Scalar &pressure,
        const std::array<Scalar, numComponents> &composition,
        const Scalar &compressibility,
        const Scalar &temperature) const
    {
        return viscosityFromMolarDensity(
            composition,
            molarDensity(pressure, compressibility, temperature),
            temperature);
    }

private:
    void initializeLbcConstants_()
    {
        constexpr double molarMassScale = 1000.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const double mwSqrt =
                std::sqrt(molarMassScale * mixture_.molecularWeight(i));
            const double tcRankine =
                mixture_.criticalTemperature(i) / units::kelvinPerRankine;
            const double pcPsi = mixture_.criticalPressure(i) / units::psi;
            lbcMolarMassSqrt_[idx] = mwSqrt;
            lbcEpsilon_[idx] =
                5.4402 * std::pow(tcRankine, 1.0 / 6.0) /
                (mwSqrt * std::pow(pcPsi, 2.0 / 3.0) * units::centipoise);
        }
    }

    CompositionalMixture<Indices> mixture_{};
    std::array<double, numComponents> lbcMolarMassSqrt_{};
    std::array<double, numComponents> lbcEpsilon_{};
};

} // namespace MPMC
