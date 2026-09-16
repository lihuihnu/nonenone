/**
 * @file aqueous_volume.hpp
 * @brief 基于 IAPWS-IF97 Region 1/2/3 与 Garcia CO2 数据的富水相摩尔体积。
 */
#pragma once

#include <natural/numerics.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace MPMC
{

/** @brief 计算 IAPWS-IF97 Region 1 压缩液态水密度。 */
struct IapwsIf97Region1 final
{
    template <class Scalar>
    [[nodiscard]] static Scalar density(
        const Scalar &pressure,
        const Scalar &temperature)
    {
        if (!(scalarValue(pressure) > 0.0) ||
            !(scalarValue(temperature) >= 273.15) ||
            !(scalarValue(temperature) <= 623.15))
        {
            throw std::invalid_argument(
                "IAPWS-IF97 Region 1 requires p>0 and 273.15<=T<=623.15 K.");
        }

        constexpr std::array<int, 34> I{
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
            2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};
        constexpr std::array<int, 34> J{
            -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1,
            3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38,
            -39, -40, -41};
        constexpr std::array<double, 34> n{
            0.14632971213167, -0.84548187169114, -3.7563603672040,
            3.3855169168385, -0.95791963387872, 0.15772038513228,
            -0.016616417199501, 0.00081214629983568, 0.00028319080123804,
            -0.00060706301565874, -0.018990068218419, -0.032529748770505,
            -0.021841717175414, -0.000052838357969930,
            -0.00047184321073267, -0.00030001780793026,
            0.000047661393906987, -0.0000044141845330846,
            -7.2694996297594e-16, -0.000031679644845054,
            -0.0000028270797985312, -8.5205128120103e-10,
            -0.0000022425281908000, -6.5171222895601e-7,
            -1.4341729937924e-13, -4.0516996860117e-7,
            -1.2734301741641e-9, -1.7424871230634e-10,
            -6.8762131295531e-19, 1.4478307828521e-20,
            2.6335781662795e-23, -1.1947622640071e-23,
            1.8228094581404e-24, -9.3537087292458e-26};

        constexpr double pressureStar = 16.53e6;
        constexpr double temperatureStar = 1386.0;
        constexpr double specificGasConstant = 461.526; // J/(kg K)
        const Scalar pi = pressure / pressureStar;
        const Scalar tau = temperatureStar / temperature;
        Scalar gammaPi = 0.0;
        for (std::size_t term = 0; term < n.size(); ++term)
        {
            if (I[term] == 0)
                continue;
            gammaPi -= n[term] * static_cast<double>(I[term]) *
                pow(Scalar(7.1) - pi, I[term] - 1) *
                pow(tau - Scalar(1.222), J[term]);
        }
        const Scalar specificVolume =
            pi * gammaPi * specificGasConstant * temperature / pressure;
        if (!(scalarValue(specificVolume) > 0.0) ||
            !std::isfinite(scalarValue(specificVolume)))
        {
            throw std::runtime_error(
                "IAPWS-IF97 Region-1 evaluation produced an invalid volume.");
        }
        return Scalar(1.0) / specificVolume;
    }
};

/** @brief 计算 IAPWS-IF97 Region 2 蒸汽/低密度超临界水密度。 */
struct IapwsIf97Region2 final
{
    template <class Scalar>
    [[nodiscard]] static Scalar density(
        const Scalar &pressure,
        const Scalar &temperature)
    {
        if (!(scalarValue(pressure) > 0.0) ||
            !(scalarValue(pressure) <= 100.0e6) ||
            !(scalarValue(temperature) >= 273.15) ||
            !(scalarValue(temperature) <= 1073.15))
        {
            throw std::invalid_argument(
                "IAPWS-IF97 Region 2 requires 0<p<=100 MPa and "
                "273.15<=T<=1073.15 K.");
        }

        constexpr std::array<int, 43> I{
            1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4,
            5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20,
            20, 20, 21, 22, 23, 24, 24, 24};
        constexpr std::array<int, 43> J{
            0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3,
            7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57,
            20, 35, 48, 21, 53, 39, 26, 40, 58};
        constexpr std::array<double, 43> n{
            -0.17731742473213e-2, -0.17834862292358e-1,
            -0.45996013696365e-1, -0.57581259083432e-1,
            -0.50325278727930e-1, -0.33032641670203e-4,
            -0.18948987516315e-3, -0.39392777243355e-2,
            -0.43797295650573e-1, -0.26674547914087e-4,
             0.20481737692309e-7,  0.43870667284435e-6,
            -0.32277677238570e-4, -0.15033924542148e-2,
            -0.40668253562649e-1, -0.78847309559367e-9,
             0.12790717852285e-7,  0.48225372718507e-6,
             0.22922076337661e-5, -0.16714766451061e-10,
            -0.21171472321355e-2, -0.23895741934104e2,
            -0.59059564324270e-17, -0.12621808899101e-5,
            -0.38946842435739e-1,  0.11256211360459e-10,
            -0.82311340897998e1,   0.19809712802088e-7,
             0.10406965210174e-18, -0.10234747095929e-12,
            -0.10018179379511e-8, -0.80882908646985e-10,
             0.10693031879409,    -0.33662250574171,
             0.89185845355421e-24, 0.30629316876232e-12,
            -0.42002467698208e-5, -0.59056029685639e-25,
             0.37826947613457e-5, -0.12768608934681e-14,
             0.73087610595061e-28, 0.55414715350778e-16,
            -0.94369707241210e-6};

        constexpr double pressureStar = 1.0e6;
        constexpr double temperatureStar = 540.0;
        constexpr double specificGasConstant = 461.526; // J/(kg K)
        const Scalar pi = pressure / pressureStar;
        const Scalar tau = temperatureStar / temperature;
        Scalar gammaPi = Scalar(1.0) / pi;
        for (std::size_t term = 0; term < n.size(); ++term)
        {
            gammaPi += n[term] * static_cast<double>(I[term]) *
                pow(pi, I[term] - 1) *
                pow(tau - Scalar(0.5), J[term]);
        }
        const Scalar specificVolume =
            pi * gammaPi * specificGasConstant * temperature / pressure;
        if (!(scalarValue(specificVolume) > 0.0) ||
            !std::isfinite(scalarValue(specificVolume)))
        {
            throw std::runtime_error(
                "IAPWS-IF97 Region-2 evaluation produced an invalid volume.");
        }
        return Scalar(1.0) / specificVolume;
    }
};

/** @brief IAPWS-IF97 Region 3 高密度/超临界水的 Helmholtz 方程与密度反算。 */
struct IapwsIf97Region3 final
{
    template <class Scalar>
    [[nodiscard]] static Scalar pressureFromDensity(
        const Scalar &density,
        const Scalar &temperature)
    {
        const auto derivatives = deltaDerivatives_(density, temperature);
        constexpr double specificGasConstant = 461.526; // J/(kg K)
        return density * specificGasConstant * temperature
            * derivatives.delta * derivatives.phiDelta;
    }

    template <class Scalar>
    [[nodiscard]] static Scalar density(
        const Scalar &pressure,
        const Scalar &temperature)
    {
        const double p = scalarValue(pressure);
        const double t = scalarValue(temperature);
        if (!(p > 0.0) || !(p <= 100.0e6) ||
            !(t >= 623.15) || !(t <= 863.15))
        {
            throw std::invalid_argument(
                "IAPWS-IF97 Region 3 requires 0<p<=100 MPa and "
                "623.15<=T<=863.15 K.");
        }

        // The physical Region-3 branch is close to 500 kg/m3 over the
        // near-critical states that trigger this closure.  A safeguarded
        // Newton solve avoids the nonphysical very-high-density polynomial
        // root while retaining the official Helmholtz equation itself.
        double rho = 500.0;
        bool converged = false;
        for (int iteration = 0; iteration < 80; ++iteration)
        {
            const double value = pressureFromDensity(rho, t) - p;
            const double derivative = pressureDerivative_(rho, t);
            if (std::isfinite(value) && std::isfinite(derivative) && derivative > 0.0)
            {
                const double candidate = rho - value / derivative;
                rho = std::clamp(candidate, 20.0, 900.0);
            }
            else
            {
                rho = 0.5 * (rho + 500.0);
            }
            if (std::abs(value) <= 1.0e-10 * std::max(p, 1.0))
            {
                converged = true;
                break;
            }
        }
        if (!converged || !(rho > 0.0) || !std::isfinite(rho))
            throw std::runtime_error("IAPWS-IF97 Region-3 density solve did not converge.");

        // Reapply Newton with Scalar arithmetic.  Starting from the converged
        // value makes the scalar residual negligible while implicit pressure
        // and temperature derivatives propagate through the same equation.
        Scalar result = rho;
        for (int iteration = 0; iteration < 3; ++iteration)
        {
            const Scalar residual =
                pressureFromDensity(result, temperature) - pressure;
            result -= residual / pressureDerivative_(result, temperature);
        }
        if (!(scalarValue(result) > 0.0) || !std::isfinite(scalarValue(result)))
            throw std::runtime_error("IAPWS-IF97 Region-3 density is invalid.");
        return result;
    }

private:
    template <class Scalar>
    struct DeltaDerivatives
    {
        Scalar delta{};
        Scalar phiDelta{};
        Scalar phiDeltaDelta{};
    };

    template <class Scalar>
    [[nodiscard]] static DeltaDerivatives<Scalar> deltaDerivatives_(
        const Scalar &density,
        const Scalar &temperature)
    {
        constexpr std::array<int, 39> I{
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3,
            3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11};
        constexpr std::array<int, 39> J{
            0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0,
            2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26,
            0, 1, 26};
        constexpr std::array<double, 39> n{
            -0.15732845290239e2,  0.20944396974307e2,
            -0.76867707878716e1,  0.26185947787954e1,
            -0.28080781148620e1,  0.12053369696517e1,
            -0.84566812812502e-2, -0.12654315477714e1,
            -0.11524407806681e1,  0.88521043984318,
            -0.64207765181607,    0.38493460186671,
            -0.85214708824206,    0.48972281541877e1,
            -0.30502617256965e1,  0.39420536879154e-1,
             0.12558408424308,   -0.27999329698710,
             0.13899799569460e1, -0.20189915023570e1,
            -0.82147637173963e-2, -0.47596035734923,
             0.43984074473500e-1, -0.44476435428739,
             0.90572070719733,    0.70522450087967,
             0.10770512626332,   -0.32913623258954,
            -0.50871062041158,   -0.22175400873096e-1,
             0.94260751665092e-1, 0.16436278447961,
            -0.13503372241348e-1, -0.14834345352472e-1,
             0.57922953628084e-3, 0.32308904703711e-2,
             0.80964802996215e-4, -0.16557679795037e-3,
            -0.44923899061815e-4};

        constexpr double criticalDensity = 322.0;
        constexpr double criticalTemperature = 647.096;
        constexpr double logarithmicCoefficient = 1.0658070028513;
        const Scalar delta = density / criticalDensity;
        const Scalar tau = criticalTemperature / temperature;
        Scalar phiDelta = logarithmicCoefficient / delta;
        Scalar phiDeltaDelta =
            -logarithmicCoefficient / (delta * delta);
        for (std::size_t term = 0; term < n.size(); ++term)
        {
            if (I[term] == 0)
                continue;
            phiDelta += n[term] * static_cast<double>(I[term])
                * pow(delta, I[term] - 1) * pow(tau, J[term]);
            if (I[term] > 1)
                phiDeltaDelta += n[term] * static_cast<double>(I[term])
                    * static_cast<double>(I[term] - 1)
                    * pow(delta, I[term] - 2) * pow(tau, J[term]);
        }
        return {delta, phiDelta, phiDeltaDelta};
    }

    template <class Scalar>
    [[nodiscard]] static Scalar pressureDerivative_(
        const Scalar &density,
        const Scalar &temperature)
    {
        const auto derivatives = deltaDerivatives_(density, temperature);
        constexpr double specificGasConstant = 461.526;
        return specificGasConstant * temperature *
            (Scalar(2.0) * derivatives.delta * derivatives.phiDelta
             + derivatives.delta * derivatives.delta
                 * derivatives.phiDeltaDelta);
    }
};

/** @brief 选择当前富水闭包已实现的 IF97 单相区域。 */
struct IapwsIf97WaterDensity final
{
    [[nodiscard]] static double region23BoundaryPressure(double temperature)
    {
        return (348.05185628969 - 1.1671859879975 * temperature
                + 0.0010192970039326 * temperature * temperature) * 1.0e6;
    }

    template <class Scalar>
    [[nodiscard]] static Scalar density(
        const Scalar &pressure,
        const Scalar &temperature)
    {
        const double p = scalarValue(pressure);
        const double t = scalarValue(temperature);
        if (t <= 623.15)
            return IapwsIf97Region1::density(pressure, temperature);
        if (t <= 863.15 && p <= region23BoundaryPressure(t))
            return IapwsIf97Region2::density(pressure, temperature);
        if (t <= 863.15 && p <= 100.0e6)
            return IapwsIf97Region3::density(pressure, temperature);
        if (t <= 1073.15 && p <= 100.0e6 && t > 863.15)
            return IapwsIf97Region2::density(pressure, temperature);
        throw std::invalid_argument("IAPWS aqueous density state is outside Regions 1-3.");
    }
};

/** @brief 计算 Garcia (2001) 水中 CO2 表观摩尔体积 [m3/mol]。 */
struct Garcia2001Co2ApparentVolume final
{
    template <class Scalar>
    [[nodiscard]] static Scalar molarVolume(const Scalar &temperature)
    {
        const Scalar celsius = temperature - 273.15;
        const Scalar celsius2 = celsius * celsius;
        const Scalar cubicCentimetresPerMole =
            37.51 - 9.585e-2 * celsius + 8.740e-4 * celsius2
            - 5.044e-7 * celsius2 * celsius;
        return cubicCentimetresPerMole * 1.0e-6;
    }
};

/**
 * @brief H2O/CO2 主导富水相的摩尔体积闭包。
 *
 * The closure changes only phase volume/density. Fugacity, Z, stability and
 * flash equations remain those of the selected thermodynamic model.
 */
template <int NumComponents>
class IapwsGarciaAqueousVolume final
{
public:
    using DoubleComposition = std::array<double, NumComponents>;

    IapwsGarciaAqueousVolume() = default;

    IapwsGarciaAqueousVolume(
        int waterComponent,
        int co2Component,
        double waterMolarMass,
        double maximumUnsupportedMoleFraction = 1.0e-4)
        : waterComponent_(waterComponent),
          co2Component_(co2Component),
          waterMolarMass_(waterMolarMass),
          maximumUnsupportedMoleFraction_(maximumUnsupportedMoleFraction)
    {
        if (waterComponent_ < 0 || waterComponent_ >= NumComponents ||
            co2Component_ < 0 || co2Component_ >= NumComponents ||
            waterComponent_ == co2Component_)
            throw std::invalid_argument("Invalid aqueous H2O/CO2 component indices.");
        if (!(waterMolarMass_ > 0.0) ||
            !(maximumUnsupportedMoleFraction_ >= 0.0))
            throw std::invalid_argument("Invalid IAPWS-Garcia configuration.");
    }

    /** @brief 判断组成是否位于 H2O+CO2 富水闭包的显式适用域内。 */
    [[nodiscard]] bool compositionSupported(
        const DoubleComposition &composition) const noexcept
    {
        double unsupported = 0.0;
        for (int component = 0; component < NumComponents; ++component)
        {
            if (component != waterComponent_ && component != co2Component_)
                unsupported += composition[static_cast<std::size_t>(component)];
        }
        return std::isfinite(unsupported) &&
            unsupported <= maximumUnsupportedMoleFraction_;
    }

    template <class Scalar>
    [[nodiscard]] Scalar molarVolume(
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
            std::ostringstream message;
            message << "IAPWS-Garcia aqueous volume cannot extrapolate a non-trace third solute"
                    << " (unsupported mole fraction=" << std::setprecision(16)
                    << scalarValue(unsupported) << ", maximum="
                    << maximumUnsupportedMoleFraction_ << ").";
            throw std::runtime_error(message.str());
        }

        const auto water = static_cast<std::size_t>(waterComponent_);
        const auto co2 = static_cast<std::size_t>(co2Component_);
        const Scalar pureWaterDensity =
            IapwsIf97WaterDensity::density(pressure, temperature);
        const Scalar volume = composition[water] * waterMolarMass_ /
                pureWaterDensity
            + composition[co2] *
                Garcia2001Co2ApparentVolume::molarVolume(temperature);
        if (!(scalarValue(volume) > 0.0) ||
            !std::isfinite(scalarValue(volume)))
            throw std::runtime_error("IAPWS-Garcia aqueous molar volume is invalid.");
        return volume;
    }

private:
    int waterComponent_{-1};
    int co2Component_{-1};
    double waterMolarMass_{0.0};
    double maximumUnsupportedMoleFraction_{1.0e-4};
};

} // namespace MPMC
