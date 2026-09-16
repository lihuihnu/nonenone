/**
 * @file aqueous_viscosity.hpp
 * @brief IAPWS 纯水与 McBride-Wright H2O-CO2 水相黏度闭包。
 */
#pragma once

#include <natural/numerics.hpp>
#include <natural/thermo/aqueous_volume.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief IAPWS-2008 工业纯水黏度公式及显式组成适用域。
 *
 * The dilute-gas and finite-density terms from IAPWS R12-08 are evaluated
 * with the production IF97 water density.  The critical-enhancement factor is
 * deliberately set to one, matching the industrial formulation outside the
 * narrow critical-enhancement region.  Dissolved non-water material is only
 * admitted up to an explicit mole-fraction limit; no unvalidated mixture
 * correction is hidden in this pure-water reference closure.
 */
template <int NumComponents>
class Iapws2008IndustrialAqueousViscosity final
{
public:
    Iapws2008IndustrialAqueousViscosity() = default;

    explicit Iapws2008IndustrialAqueousViscosity(
        int waterComponent,
        double maximumSoluteMoleFraction = 0.02)
        : waterComponent_(waterComponent),
          maximumSoluteMoleFraction_(maximumSoluteMoleFraction)
    {
        if (waterComponent_ < 0 || waterComponent_ >= NumComponents ||
            !(maximumSoluteMoleFraction_ >= 0.0))
            throw std::invalid_argument(
                "Invalid IAPWS-2008 aqueous-viscosity configuration.");
    }

    template <class Scalar>
    [[nodiscard]] Scalar viscosity(
        const Scalar &pressure,
        const Scalar &temperature,
        const std::array<Scalar, NumComponents> &composition) const
    {
        Scalar solute = 0.0;
        for (int component = 0; component < NumComponents; ++component)
            if (component != waterComponent_)
                solute += composition[static_cast<std::size_t>(component)];
        if (scalarValue(solute) > maximumSoluteMoleFraction_)
            throw std::runtime_error(
                "IAPWS-2008 pure-water viscosity cannot extrapolate the active water phase's solute fraction.");

        const Scalar density =
            IapwsIf97WaterDensity::density(pressure, temperature);
        const Scalar reducedTemperature = temperature / 647.096;
        const Scalar reducedDensity = density / 322.0;
        constexpr std::array<double, 4> h{
            1.67752, 2.20462, 0.6366564, -0.241605};
        Scalar denominator = 0.0;
        Scalar temperaturePower = 1.0;
        for (std::size_t i = 0; i < h.size(); ++i)
        {
            if (i > 0)
                temperaturePower *= reducedTemperature;
            denominator += h[i] / temperaturePower;
        }
        const Scalar mu0 = 100.0 * sqrt(reducedTemperature) / denominator;

        struct Coefficient
        {
            int i;
            int j;
            double value;
        };
        constexpr std::array<Coefficient, 21> coefficients{{
            {0, 0, 5.20094e-1}, {1, 0, 8.50895e-2}, {2, 0, -1.08374},
            {3, 0, -2.89555e-1}, {0, 1, 2.22531e-1}, {1, 1, 9.99115e-1},
            {2, 1, 1.88797}, {3, 1, 1.26613}, {5, 1, 1.20573e-1},
            {0, 2, -2.81378e-1}, {1, 2, -9.06851e-1}, {2, 2, -7.72479e-1},
            {3, 2, -4.89837e-1}, {4, 2, -2.57040e-1}, {0, 3, 1.61913e-1},
            {1, 3, 2.57399e-1}, {0, 4, -3.25372e-2}, {3, 4, 6.98452e-2},
            {4, 5, 8.72102e-3}, {3, 6, -4.35673e-3}, {5, 6, -5.93264e-4}}};
        const Scalar temperatureVariable = 1.0 / reducedTemperature - 1.0;
        const Scalar densityVariable = reducedDensity - 1.0;
        Scalar exponentSeries = 0.0;
        for (const auto &coefficient : coefficients)
            exponentSeries += coefficient.value
                * pow(temperatureVariable, coefficient.i)
                * pow(densityVariable, coefficient.j);
        const Scalar mu1 = exp(reducedDensity * exponentSeries);
        const Scalar result = mu0 * mu1 * 1.0e-6;
        if (!(scalarValue(result) > 0.0) ||
            !std::isfinite(scalarValue(result)))
            throw std::runtime_error("IAPWS-2008 viscosity is invalid.");
        return result;
    }

private:
    int waterComponent_{-1};
    double maximumSoluteMoleFraction_{0.02};
};

template <int NumComponents>
class McBrideWright2015AqueousViscosity final
{
public:
    McBrideWright2015AqueousViscosity() = default;

    McBrideWright2015AqueousViscosity(
        int waterComponent,
        int co2Component,
        double maximumUnsupportedMoleFraction = 1.0e-4)
        : waterComponent_(waterComponent),
          co2Component_(co2Component),
          maximumUnsupportedMoleFraction_(maximumUnsupportedMoleFraction)
    {
        if (waterComponent_ < 0 || waterComponent_ >= NumComponents ||
            co2Component_ < 0 || co2Component_ >= NumComponents ||
            waterComponent_ == co2Component_ ||
            !(maximumUnsupportedMoleFraction_ >= 0.0))
        {
            throw std::invalid_argument(
                "Invalid McBride-Wright aqueous-viscosity configuration.");
        }
    }

    template <class Scalar>
    [[nodiscard]] Scalar viscosity(
        const Scalar &pressure,
        const Scalar &temperature,
        const std::array<Scalar, NumComponents> &composition) const
    {
        Scalar unsupported = 0.0;
        for (int component = 0; component < NumComponents; ++component)
        {
            if (component != waterComponent_ && component != co2Component_)
                unsupported += composition[static_cast<std::size_t>(component)];
        }
        if (scalarValue(unsupported) > maximumUnsupportedMoleFraction_)
        {
            throw std::runtime_error(
                "McBride-Wright viscosity cannot extrapolate a non-trace third solute.");
        }

        constexpr double a = -3.705013;
        constexpr double b = 0.00289258;
        constexpr double c = 3.98950;
        constexpr double d = -0.00326;
        constexpr double e1 = 65.55968;
        constexpr double e2 = 2.46811;
        constexpr double t0 = 141.5; // K
        constexpr double p0 = 1.0e6; // Pa
        const Scalar reducedPressure = pressure / p0;
        const Scalar shiftedTemperature = temperature / t0 - 1.0;
        const Scalar xCo2 = composition[static_cast<std::size_t>(co2Component_)];
        const Scalar logMilliPascalSeconds =
            a + b * reducedPressure +
            (c + d * reducedPressure) / shiftedTemperature +
            e1 * exp(-e2 * shiftedTemperature) * xCo2;
        const Scalar result = exp(logMilliPascalSeconds) * 1.0e-3;
        if (!(scalarValue(result) > 0.0) || !std::isfinite(scalarValue(result)))
            throw std::runtime_error("McBride-Wright viscosity is invalid.");
        return result; // Pa s
    }

private:
    int waterComponent_{-1};
    int co2Component_{-1};
    double maximumUnsupportedMoleFraction_{1.0e-4};
};

} // namespace MPMC
