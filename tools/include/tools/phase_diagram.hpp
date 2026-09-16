/**
 * @file phase_diagram.hpp
 * @brief 基于 production flash 的 P–T、组成和三元相图扫描工具。
 */
#pragma once

#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace MPMC::tools
{

/**
 * @brief 普通水临界点的 IAPWS 参考坐标。
 *
 * These values are plotting/reference metadata only; no IAPWS water-property
 * equation is substituted into PR, Soreide-Whitson, or CPA calculations.
 */
struct WaterCriticalReference final
{
    static constexpr double temperatureK = 647.096;
    static constexpr double pressurePa = 22.064e6;
    static constexpr double pressureBar = pressurePa / 1.0e5;
};

/** @brief 一维扫描轴使用的取点间距类型。 */
enum class AxisSpacing
{
    Linear,
    Logarithmic
};

/**
 * @brief 包含两端点的一维扫描轴。
 *
 * `points` includes both endpoints.  Logarithmic spacing requires positive
 * endpoints and is useful for wide pressure ranges.
 */
struct ScanAxis
{
    double minimum{0.0};
    double maximum{0.0};
    std::size_t points{1};
    AxisSpacing spacing{AxisSpacing::Linear};

    [[nodiscard]] std::vector<double> values() const
    {
        if (points == 0)
            throw std::invalid_argument("Phase-diagram scan axis needs at least one point.");
        if (!std::isfinite(minimum) || !std::isfinite(maximum))
            throw std::invalid_argument("Phase-diagram scan axis endpoints must be finite.");
        if (maximum < minimum)
            throw std::invalid_argument("Phase-diagram scan axis maximum is smaller than minimum.");
        if (spacing == AxisSpacing::Logarithmic && (!(minimum > 0.0) || !(maximum > 0.0)))
            throw std::invalid_argument("Logarithmic phase-diagram axis endpoints must be positive.");

        std::vector<double> result(points, minimum);
        if (points == 1)
            return result;

        if (spacing == AxisSpacing::Linear)
        {
            const double step = (maximum - minimum) / static_cast<double>(points - 1);
            for (std::size_t i = 0; i < points; ++i)
                result[i] = minimum + step * static_cast<double>(i);
        }
        else
        {
            const double logMin = std::log(minimum);
            const double logMax = std::log(maximum);
            const double step = (logMax - logMin) / static_cast<double>(points - 1);
            for (std::size_t i = 0; i < points; ++i)
                result[i] = std::exp(logMin + step * static_cast<double>(i));
        }
        return result;
    }
};

/** @brief 指定扫描允许生成全部相还是限制在给定 active set。 */
enum class FlashMode
{
    Unrestricted,
    Restricted
};

/**
 * @brief 相图采样使用的 flash 策略。
 *
 * Restricted mode is useful for conventional hydrocarbon VLE diagrams where a
 * second PR liquid slot must not be interpreted as a distinct water-rich phase.
 * Fully compositional H2O systems should normally use Unrestricted mode.
 */
struct PhaseDiagramFlashPolicy
{
    FlashMode mode{FlashMode::Unrestricted};
    PhasePresence allowedPhases{PhasePresence::all()};
    bool continueAfterFailure{true};

    [[nodiscard]] static PhaseDiagramFlashPolicy unrestricted() noexcept
    {
        return {};
    }

    [[nodiscard]] static PhaseDiagramFlashPolicy restricted(PhasePresence allowed) noexcept
    {
        PhaseDiagramFlashPolicy policy;
        policy.mode = FlashMode::Restricted;
        policy.allowedPhases = allowed;
        return policy;
    }
};

/**
 * @brief 相包络提取的边界识别配置。
 *
 * The default interpretation matches conventional hydrocarbon vapor-liquid
 * diagrams: single liquid (code 1), single vapor (code 2), and two-phase O+G
 * (code 3).  For more specialised studies these codes may be reassigned to a
 * different pair of single/two-phase states as long as the one-dimensional
 * pressure slices remain monotone.
 */
struct EnvelopeOptions
{
    int liquidPhaseCode{1};
    int vaporPhaseCode{2};
    int twoPhaseCode{3};
    int maxRefinementIterations{40};
    double relativePressureTolerance{1.0e-8};
    // Quality gate for the envelope-coalescence estimate.  The score is
    // relative bubble/dew pressure span + 0.5*L1(oil,gas composition).
    double criticalEstimateScoreTolerance{0.20};
};

/** @brief 一个 P-T-z flash 样本及其扫描元数据。 */
template <class Indices>
struct PhaseDiagramSample
{
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;
    using FlashResult = ThreePhaseFlashResult<Indices>;

    double pressurePa{0.0};
    double temperatureK{0.0};
    Composition overallComposition{};
    FlashResult flash{};
    int statusCode{0}; // 0=converged, 1=non-converged, 2=exception

    [[nodiscard]] int phaseCode() const noexcept
    {
        return statusCode == 0 && flash.converged
            ? static_cast<int>(flash.presence.bits())
            : 0;
    }

    [[nodiscard]] int phaseCount() const noexcept
    {
        return phaseCode() == 0 ? 0 : flash.presence.count();
    }
};

/** @brief 规则压力-温度相图数据。 */
template <class Indices>
struct PressureTemperatureMap
{
    using Composition = std::array<double, Indices::numComponents>;
    ScanAxis temperatureAxis{};
    ScanAxis pressureAxis{};
    Composition overallComposition{};
    std::vector<PhaseDiagramSample<Indices>> samples{};
};

/** @brief 沿两个总体组成之间线性路径的压力扫描数据。 */
template <class Indices>
struct PressureCompositionMap
{
    using Composition = std::array<double, Indices::numComponents>;
    double temperatureK{0.0};
    ScanAxis pressureAxis{};
    ScanAxis pathAxis{0.0, 1.0, 2, AxisSpacing::Linear};
    Composition startComposition{};
    Composition endComposition{};
    std::vector<double> pathFraction{};
    std::vector<PhaseDiagramSample<Indices>> samples{};
};

/** @brief 固定 P、T 的三元组成相图。 */
template <class Indices>
struct TernaryCompositionMap
{
    double pressurePa{0.0};
    double temperatureK{0.0};
    std::array<int, 3> components{0, 1, 2};
    std::size_t subdivisions{1};
    std::vector<std::array<double, 3>> ternaryFraction{};
    std::vector<PhaseDiagramSample<Indices>> samples{};
};

/** @brief 泡点/露点包络线上的一个坐标点。 */
struct EnvelopePoint
{
    double coordinate{0.0};
    bool hasDew{false};
    double dewPressurePa{std::numeric_limits<double>::quiet_NaN()};
    bool hasBubble{false};
    double bubblePressurePa{std::numeric_limits<double>::quiet_NaN()};
    int boundedTwoPhaseIntervals{0};

    [[nodiscard]] double envelopeWidthPa() const noexcept
    {
        return (hasDew && hasBubble) ? bubblePressurePa - dewPressurePa
                                     : std::numeric_limits<double>::quiet_NaN();
    }
};

/** @brief 固定总体组成的泡点/露点压力曲线。 */
template <class Indices>
struct PressureTemperatureEnvelope
{
    using Composition = std::array<double, Indices::numComponents>;
    ScanAxis temperatureAxis{};
    ScanAxis pressureAxis{};
    Composition overallComposition{};
    EnvelopeOptions options{};
    std::vector<EnvelopePoint> points{};
    bool hasCriticalPoint{false};
    double criticalTemperatureK{std::numeric_limits<double>::quiet_NaN()};
    double criticalPressurePa{std::numeric_limits<double>::quiet_NaN()};
    double criticalEnvelopeWidthPa{std::numeric_limits<double>::quiet_NaN()};
    double criticalRelativeEnvelopeWidth{std::numeric_limits<double>::quiet_NaN()};
    double criticalOilGasCompositionL1{std::numeric_limits<double>::quiet_NaN()};
    double criticalEstimateScore{std::numeric_limits<double>::quiet_NaN()};
    bool criticalEstimateWithinTolerance{false};
};

/** @brief 定温组成路径上的泡点/露点压力曲线。 */
template <class Indices>
struct PressureCompositionEnvelope
{
    using Composition = std::array<double, Indices::numComponents>;
    double temperatureK{0.0};
    ScanAxis pressureAxis{};
    ScanAxis pathAxis{};
    Composition startComposition{};
    Composition endComposition{};
    EnvelopeOptions options{};
    std::vector<EnvelopePoint> points{};
};

/** @brief 组成路径临界轨迹上的一个近似临界点。 */
template <class Indices>
struct CriticalLocusPoint
{
    using Composition = std::array<double, Indices::numComponents>;
    double pathFraction{0.0};
    Composition overallComposition{};
    bool valid{false};
    double criticalTemperatureK{std::numeric_limits<double>::quiet_NaN()};
    double criticalPressurePa{std::numeric_limits<double>::quiet_NaN()};
    double criticalEnvelopeWidthPa{std::numeric_limits<double>::quiet_NaN()};
    double criticalRelativeEnvelopeWidth{std::numeric_limits<double>::quiet_NaN()};
    double criticalOilGasCompositionL1{std::numeric_limits<double>::quiet_NaN()};
    double criticalEstimateScore{std::numeric_limits<double>::quiet_NaN()};
    bool criticalEstimateWithinTolerance{false};
};

/**
 * @brief 沿组成路径得到的近似临界轨迹。
 *
 * Each point is obtained by first extracting a fixed-composition P-T envelope
 * on the provided scan grid and then selecting the smallest bubble-dew span as
 * an approximate critical state.  This is intended for paper-quality trend and
 * teaching plots; it is not a full critical continuation solver.
 */
template <class Indices>
struct CriticalLocusPath
{
    using Composition = std::array<double, Indices::numComponents>;
    ScanAxis pathAxis{};
    Composition startComposition{};
    Composition endComposition{};
    ScanAxis temperatureAxis{};
    ScanAxis pressureAxis{};
    EnvelopeOptions options{};
    std::vector<CriticalLocusPoint<Indices>> points{};
};

/** @brief 非限制 O/G/W P-T 图中的一次相存在位切换。 */
enum class MultiphaseBoundaryKind
{
    OilOnset = 1,
    GasOnset = 2,
    WaterOnset = 3
};

/** @brief O/G/W 相出现边界细化的数值控制参数。 */
struct MultiphaseBoundaryOptions
{
    int maxRefinementIterations{32};
    double relativePressureTolerance{1.0e-7};
};

template <class Indices>
struct MultiphaseBoundaryPoint
{
    double temperatureK{0.0};
    double pressurePa{0.0};
    MultiphaseBoundaryKind kind{MultiphaseBoundaryKind::WaterOnset};
    int lowerPhaseCode{0};
    int upperPhaseCode{0};
};

template <class Indices>
struct PressureTemperaturePhaseBoundaries
{
    using Composition = std::array<double, Indices::numComponents>;
    ScanAxis temperatureAxis{};
    ScanAxis pressureAxis{};
    Composition overallComposition{};
    std::vector<MultiphaseBoundaryPoint<Indices>> points{};
};

namespace detail
{

template <std::size_t N>
[[nodiscard]] inline std::array<double, N> normalizedComposition(std::array<double, N> z)
{
    double sum = 0.0;
    for (double value : z)
    {
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument("Phase-diagram overall composition must be finite and non-negative.");
        sum += value;
    }
    if (!(sum > 0.0))
        throw std::invalid_argument("Phase-diagram overall composition must have positive sum.");
    for (double &value : z)
        value /= sum;
    return z;
}

[[nodiscard]] inline std::string sanitizeColumnName(std::string_view name)
{
    std::string result;
    result.reserve(name.size());
    for (const char ch : name)
    {
        const bool alphaNumeric =
            (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
            (ch >= '0' && ch <= '9');
        result.push_back(alphaNumeric ? ch : '_');
    }
    if (result.empty())
        result = "component";
    return result;
}

template <class Indices, class Name>
inline void writeCommonHeader(
    std::ostream &out,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    out << "pressure_Pa,pressure_bar,temperature_K,status_code,phase_code,phase_count,iterations"
        << ",beta_o,beta_g,beta_w,So,Sg,Sw,Z_o,Z_g,Z_w,rho_m_o,rho_m_g,rho_m_w";
    for (const auto &name : componentNames)
        out << ",z_" << sanitizeColumnName(name);
    for (const char *prefix : {"xo_", "yg_", "xw_"})
        for (const auto &name : componentNames)
            out << ',' << prefix << sanitizeColumnName(name);
    for (const char *prefix : {"Kg_", "Kw_"})
        for (const auto &name : componentNames)
            out << ',' << prefix << sanitizeColumnName(name);
}

template <class Indices>
inline void writeCommonRow(std::ostream &out, const PhaseDiagramSample<Indices> &sample)
{
    const auto &r = sample.flash;
    out << sample.pressurePa << ',' << sample.pressurePa / 1.0e5 << ',' << sample.temperatureK
        << ',' << sample.statusCode << ',' << sample.phaseCode() << ',' << sample.phaseCount()
        << ',' << r.iterations;
    for (double value : r.phaseMoleFraction)
        out << ',' << value;
    for (double value : r.saturation)
        out << ',' << value;
    for (double value : r.compressibility)
        out << ',' << value;
    for (double value : r.molarDensity)
        out << ',' << value;
    for (double value : sample.overallComposition)
        out << ',' << value;
    for (const auto &phaseComposition : r.composition)
        for (double value : phaseComposition)
            out << ',' << value;
    for (double value : r.vaporOilK)
        out << ',' << value;
    for (double value : r.waterOilK)
        out << ',' << value;
}

inline void prepareOutputFile(const std::filesystem::path &path)
{
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
}

[[nodiscard]] inline bool isFinitePositive(double value) noexcept
{
    return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] inline bool phaseMatches(const int code, const int expected) noexcept
{
    return code == expected;
}

[[nodiscard]] inline bool intervalSmallEnough(double a, double b, double relTol) noexcept
{
    const double width = std::abs(b - a);
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return width <= relTol * scale;
}

template <std::size_t N>
[[nodiscard]] inline std::array<double, N> interpolatedComposition(
    const std::array<double, N> &start,
    const std::array<double, N> &end,
    double fraction)
{
    std::array<double, N> z{};
    for (std::size_t i = 0; i < N; ++i)
        z[i] = (1.0 - fraction) * start[i] + fraction * end[i];
    return normalizedComposition(z);
}

template <class Indices>
[[nodiscard]] inline double oilGasCompositionL1(const PhaseDiagramSample<Indices> &sample)
{
    const int code = sample.phaseCode();
    if ((code & static_cast<int>(PhasePresence::oilBit)) == 0 ||
        (code & static_cast<int>(PhasePresence::gasBit)) == 0)
        return std::numeric_limits<double>::infinity();
    double distance = 0.0;
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const auto c = static_cast<std::size_t>(component);
        distance += std::abs(sample.flash.composition[0][c] - sample.flash.composition[1][c]);
    }
    return distance;
}

[[nodiscard]] inline MultiphaseBoundaryKind boundaryKindFromChangedBit(std::uint8_t bit)
{
    if (bit == PhasePresence::oilBit)
        return MultiphaseBoundaryKind::OilOnset;
    if (bit == PhasePresence::gasBit)
        return MultiphaseBoundaryKind::GasOnset;
    return MultiphaseBoundaryKind::WaterOnset;
}

} // namespace detail

/**
 * @brief 基于 production 三相 flash 构建的可复用相图采样器。
 *
 * The sampler never implements its own EOS or phase-equilibrium equations.  It
 * only varies P/T/z and stores the production flash result, so generated plots
 * use the same thermodynamic kernel as the reservoir simulator.
 */
template <class Indices>
class PhaseDiagramSampler final
{
public:
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;
    using Flash = CubicThreePhaseFlash<Indices>;
    using Sample = PhaseDiagramSample<Indices>;

    explicit PhaseDiagramSampler(
        const Flash &flash,
        PhaseDiagramFlashPolicy policy = PhaseDiagramFlashPolicy::unrestricted())
        : flash_(flash), policy_(policy)
    {
        if (policy_.mode == FlashMode::Restricted && policy_.allowedPhases.empty())
            throw std::invalid_argument("Restricted phase-diagram flash needs at least one allowed phase.");
    }

    [[nodiscard]] Sample evaluate(double pressurePa, double temperatureK, Composition z) const
    {
        if (!(pressurePa > 0.0) || !(temperatureK > 0.0) ||
            !std::isfinite(pressurePa) || !std::isfinite(temperatureK))
        {
            throw std::invalid_argument("Phase-diagram P and T must be finite and positive.");
        }

        Sample sample;
        sample.pressurePa = pressurePa;
        sample.temperatureK = temperatureK;
        sample.overallComposition = detail::normalizedComposition(z);

        try
        {
            sample.flash = policy_.mode == FlashMode::Restricted
                ? flash_.flashRestricted(
                    pressurePa, temperatureK, sample.overallComposition, policy_.allowedPhases)
                : flash_.flash(pressurePa, temperatureK, sample.overallComposition);
            sample.statusCode = sample.flash.converged ? 0 : 1;
            if (sample.statusCode != 0 && !policy_.continueAfterFailure)
                throw std::runtime_error("Phase-diagram flash did not converge.");
        }
        catch (const std::exception &)
        {
            if (!policy_.continueAfterFailure)
                throw;
            sample.statusCode = 2;
            sample.flash.converged = false;
            sample.flash.presence = PhasePresence(static_cast<std::uint8_t>(0));
        }
        return sample;
    }

    [[nodiscard]] PressureTemperatureMap<Indices> pressureTemperature(
        Composition z,
        ScanAxis temperatureAxis,
        ScanAxis pressureAxis) const
    {
        PressureTemperatureMap<Indices> map;
        map.temperatureAxis = temperatureAxis;
        map.pressureAxis = pressureAxis;
        map.overallComposition = detail::normalizedComposition(z);

        const auto temperatures = temperatureAxis.values();
        const auto pressures = pressureAxis.values();
        map.samples.reserve(temperatures.size() * pressures.size());
        for (double temperature : temperatures)
            for (double pressure : pressures)
                map.samples.push_back(evaluate(pressure, temperature, map.overallComposition));
        return map;
    }

    [[nodiscard]] PressureCompositionMap<Indices> pressureComposition(
        double temperatureK,
        ScanAxis pressureAxis,
        Composition startComposition,
        Composition endComposition,
        ScanAxis pathAxis = ScanAxis{0.0, 1.0, 51, AxisSpacing::Linear}) const
    {
        if (pathAxis.minimum < 0.0 || pathAxis.maximum > 1.0)
            throw std::invalid_argument("Pressure-composition path fraction must remain in [0,1].");

        PressureCompositionMap<Indices> map;
        map.temperatureK = temperatureK;
        map.pressureAxis = pressureAxis;
        map.pathAxis = pathAxis;
        map.startComposition = detail::normalizedComposition(startComposition);
        map.endComposition = detail::normalizedComposition(endComposition);

        const auto pressures = pressureAxis.values();
        const auto fractions = pathAxis.values();
        map.samples.reserve(pressures.size() * fractions.size());
        map.pathFraction.reserve(pressures.size() * fractions.size());
        for (double fraction : fractions)
        {
            const Composition z = detail::interpolatedComposition(map.startComposition, map.endComposition, fraction);
            for (double pressure : pressures)
            {
                map.pathFraction.push_back(fraction);
                map.samples.push_back(evaluate(pressure, temperatureK, z));
            }
        }
        return map;
    }

    [[nodiscard]] TernaryCompositionMap<Indices> ternaryComposition(
        double pressurePa,
        double temperatureK,
        std::array<int, 3> components,
        std::size_t subdivisions) const
    {
        if (subdivisions == 0)
            throw std::invalid_argument("Ternary phase diagram needs at least one subdivision.");
        for (int c : components)
            if (c < 0 || c >= N)
                throw std::invalid_argument("Ternary phase-diagram component index is out of range.");
        if (components[0] == components[1] || components[0] == components[2] ||
            components[1] == components[2])
            throw std::invalid_argument("Ternary phase-diagram components must be distinct.");

        TernaryCompositionMap<Indices> map;
        map.pressurePa = pressurePa;
        map.temperatureK = temperatureK;
        map.components = components;
        map.subdivisions = subdivisions;
        const std::size_t pointCount = (subdivisions + 1) * (subdivisions + 2) / 2;
        map.samples.reserve(pointCount);
        map.ternaryFraction.reserve(pointCount);

        for (std::size_t ia = 0; ia <= subdivisions; ++ia)
        {
            for (std::size_t ib = 0; ib + ia <= subdivisions; ++ib)
            {
                const std::size_t ic = subdivisions - ia - ib;
                const std::array<double, 3> fraction{
                    static_cast<double>(ia) / static_cast<double>(subdivisions),
                    static_cast<double>(ib) / static_cast<double>(subdivisions),
                    static_cast<double>(ic) / static_cast<double>(subdivisions)};
                Composition z{};
                z[static_cast<std::size_t>(components[0])] = fraction[0];
                z[static_cast<std::size_t>(components[1])] = fraction[1];
                z[static_cast<std::size_t>(components[2])] = fraction[2];
                map.ternaryFraction.push_back(fraction);
                map.samples.push_back(evaluate(pressurePa, temperatureK, z));
            }
        }
        return map;
    }

    [[nodiscard]] PressureTemperatureEnvelope<Indices> pressureTemperatureEnvelope(
        Composition z,
        ScanAxis temperatureAxis,
        ScanAxis pressureAxis,
        EnvelopeOptions options = EnvelopeOptions{}) const
    {
        const auto map = pressureTemperature(z, temperatureAxis, pressureAxis);
        return pressureTemperatureEnvelope(map, options);
    }

    [[nodiscard]] PressureTemperatureEnvelope<Indices> pressureTemperatureEnvelope(
        const PressureTemperatureMap<Indices> &map,
        EnvelopeOptions options = EnvelopeOptions{}) const
    {
        if (options.maxRefinementIterations <= 0 ||
            !(options.relativePressureTolerance > 0.0) ||
            !std::isfinite(options.relativePressureTolerance) ||
            !(options.criticalEstimateScoreTolerance > 0.0) ||
            !std::isfinite(options.criticalEstimateScoreTolerance))
            throw std::invalid_argument(
                "Envelope refinement and critical-estimate tolerances must be finite and positive.");

        PressureTemperatureEnvelope<Indices> envelope;
        envelope.temperatureAxis = map.temperatureAxis;
        envelope.pressureAxis = map.pressureAxis;
        envelope.overallComposition = map.overallComposition;
        envelope.options = options;

        const std::size_t nT = map.temperatureAxis.points;
        const std::size_t nP = map.pressureAxis.points;
        const auto temperatures = map.temperatureAxis.values();
        const auto pressures = map.pressureAxis.values();
        envelope.points.reserve(nT);

        double bestScore = std::numeric_limits<double>::infinity();
        for (std::size_t it = 0; it < nT; ++it)
        {
            std::vector<Sample> line;
            line.reserve(nP);
            for (std::size_t ip = 0; ip < nP; ++ip)
                line.push_back(map.samples[it * nP + ip]);

            EnvelopePoint point = extractEnvelopePoint(temperatures[it], map.overallComposition, pressures, line, options);
            envelope.points.push_back(point);
            const double width = point.envelopeWidthPa();
            if (!(std::isfinite(width) && width >= 0.0))
                continue;

            const double criticalPressure = 0.5 * (point.bubblePressurePa + point.dewPressurePa);
            const double relativeWidth = width / std::max(criticalPressure, 1.0);
            const Sample criticalSample = evaluate(criticalPressure, temperatures[it], map.overallComposition);
            const double compositionDistance = detail::oilGasCompositionL1(criticalSample);
            const double finiteCompositionDistance =
                std::isfinite(compositionDistance) ? compositionDistance : 2.0;
            const double score = relativeWidth + 0.5 * finiteCompositionDistance;
            if (score < bestScore)
            {
                bestScore = score;
                envelope.hasCriticalPoint = true;
                envelope.criticalTemperatureK = temperatures[it];
                envelope.criticalPressurePa = criticalPressure;
                envelope.criticalEnvelopeWidthPa = width;
                envelope.criticalRelativeEnvelopeWidth = relativeWidth;
                envelope.criticalOilGasCompositionL1 = compositionDistance;
                envelope.criticalEstimateScore = score;
                envelope.criticalEstimateWithinTolerance =
                    score <= options.criticalEstimateScoreTolerance;
            }
        }
        return envelope;
    }

    [[nodiscard]] PressureCompositionEnvelope<Indices> pressureCompositionEnvelope(
        double temperatureK,
        ScanAxis pressureAxis,
        Composition startComposition,
        Composition endComposition,
        ScanAxis pathAxis = ScanAxis{0.0, 1.0, 51, AxisSpacing::Linear},
        EnvelopeOptions options = EnvelopeOptions{}) const
    {
        const auto map = pressureComposition(temperatureK, pressureAxis, startComposition, endComposition, pathAxis);
        return pressureCompositionEnvelope(map, options);
    }

    [[nodiscard]] PressureCompositionEnvelope<Indices> pressureCompositionEnvelope(
        const PressureCompositionMap<Indices> &map,
        EnvelopeOptions options = EnvelopeOptions{}) const
    {
        PressureCompositionEnvelope<Indices> envelope;
        envelope.temperatureK = map.temperatureK;
        envelope.pressureAxis = map.pressureAxis;
        envelope.pathAxis = map.pathAxis;
        envelope.startComposition = map.startComposition;
        envelope.endComposition = map.endComposition;
        envelope.options = options;

        const std::size_t nX = map.pathAxis.points;
        const std::size_t nP = map.pressureAxis.points;
        const auto fractions = map.pathAxis.values();
        const auto pressures = map.pressureAxis.values();
        envelope.points.reserve(nX);

        for (std::size_t ix = 0; ix < nX; ++ix)
        {
            std::vector<Sample> line;
            line.reserve(nP);
            for (std::size_t ip = 0; ip < nP; ++ip)
                line.push_back(map.samples[ix * nP + ip]);
            const Composition z = detail::interpolatedComposition(map.startComposition, map.endComposition, fractions[ix]);
            EnvelopePoint point = extractEnvelopePoint(map.temperatureK, z, pressures, line, options);
            point.coordinate = fractions[ix];
            envelope.points.push_back(point);
        }
        return envelope;
    }

    [[nodiscard]] CriticalLocusPath<Indices> criticalLocus(
        ScanAxis pathAxis,
        Composition startComposition,
        Composition endComposition,
        ScanAxis temperatureAxis,
        ScanAxis pressureAxis,
        EnvelopeOptions options = EnvelopeOptions{}) const
    {
        if (pathAxis.minimum < 0.0 || pathAxis.maximum > 1.0)
            throw std::invalid_argument("Critical-locus composition path fraction must remain in [0,1].");

        CriticalLocusPath<Indices> locus;
        locus.pathAxis = pathAxis;
        locus.startComposition = detail::normalizedComposition(startComposition);
        locus.endComposition = detail::normalizedComposition(endComposition);
        locus.temperatureAxis = temperatureAxis;
        locus.pressureAxis = pressureAxis;
        locus.options = options;

        const auto fractions = pathAxis.values();
        locus.points.reserve(fractions.size());
        for (double fraction : fractions)
        {
            const Composition z = detail::interpolatedComposition(locus.startComposition, locus.endComposition, fraction);
            const auto envelope = pressureTemperatureEnvelope(z, temperatureAxis, pressureAxis, options);
            CriticalLocusPoint<Indices> point;
            point.pathFraction = fraction;
            point.overallComposition = z;
            point.valid = envelope.hasCriticalPoint;
            if (point.valid)
            {
                point.criticalTemperatureK = envelope.criticalTemperatureK;
                point.criticalPressurePa = envelope.criticalPressurePa;
                point.criticalEnvelopeWidthPa = envelope.criticalEnvelopeWidthPa;
                point.criticalRelativeEnvelopeWidth = envelope.criticalRelativeEnvelopeWidth;
                point.criticalOilGasCompositionL1 = envelope.criticalOilGasCompositionL1;
                point.criticalEstimateScore = envelope.criticalEstimateScore;
                point.criticalEstimateWithinTolerance = envelope.criticalEstimateWithinTolerance;
            }
            locus.points.push_back(point);
        }
        return locus;
    }

    /**
     * @brief 从非限制 P-T 图提取油/气/水相出现边界。
     *
     * Every changed O/G/W presence bit between adjacent converged samples is
     * refined independently.  This also resolves coarse-grid jumps such as
     * O -> O+G+W without silently discarding the gas- and water-onset events.
     * Failed samples are never bridged, so missing curve segments remain an
     * explicit numerical diagnostic rather than an inferred phase boundary.
     */
    [[nodiscard]] PressureTemperaturePhaseBoundaries<Indices>
    pressureTemperaturePhaseBoundaries(
        const PressureTemperatureMap<Indices> &map,
        MultiphaseBoundaryOptions options = {}) const
    {
        if (options.maxRefinementIterations <= 0 ||
            !(options.relativePressureTolerance > 0.0) ||
            !std::isfinite(options.relativePressureTolerance))
            throw std::invalid_argument(
                "Multiphase boundary refinement controls must be finite and positive.");

        PressureTemperaturePhaseBoundaries<Indices> result;
        result.temperatureAxis = map.temperatureAxis;
        result.pressureAxis = map.pressureAxis;
        result.overallComposition = map.overallComposition;
        const auto temperatures = map.temperatureAxis.values();
        const auto pressures = map.pressureAxis.values();
        const std::size_t nP = pressures.size();

        for (std::size_t it = 0; it < temperatures.size(); ++it)
        {
            for (std::size_t ip = 0; ip + 1 < nP; ++ip)
            {
                const Sample &a = map.samples[it * nP + ip];
                const Sample &b = map.samples[it * nP + ip + 1];
                if (a.statusCode != 0 || b.statusCode != 0)
                    continue;
                const int codeA = a.phaseCode();
                const int codeB = b.phaseCode();
                const auto changed = static_cast<std::uint8_t>(codeA ^ codeB);
                for (const std::uint8_t bit : {
                         PhasePresence::oilBit,
                         PhasePresence::gasBit,
                         PhasePresence::waterBit})
                {
                    if ((changed & bit) == 0u)
                        continue;
                    const double pressure = refinePresenceBitTransition(
                        map.overallComposition,
                        temperatures[it],
                        pressures[ip],
                        pressures[ip + 1],
                        bit,
                        options);
                    result.points.push_back({
                        temperatures[it],
                        pressure,
                        detail::boundaryKindFromChangedBit(bit),
                        codeA,
                        codeB});
                }
            }
        }
        return result;
    }

private:
    [[nodiscard]] EnvelopePoint extractEnvelopePoint(
        double temperatureK,
        const Composition &z,
        const std::vector<double> &pressures,
        const std::vector<Sample> &line,
        const EnvelopeOptions &options) const
    {
        EnvelopePoint point;
        point.coordinate = temperatureK;
        if (pressures.size() != line.size() || pressures.size() < 2)
            return point;

        for (std::size_t i = 0; i + 1 < line.size(); ++i)
        {
            const Sample &a = line[i];
            const Sample &b = line[i + 1];
            if (a.statusCode != 0 || b.statusCode != 0)
                continue;

            const int codeA = a.phaseCode();
            const int codeB = b.phaseCode();
            if (detail::phaseMatches(codeA, options.vaporPhaseCode) && detail::phaseMatches(codeB, options.twoPhaseCode))
            {
                point.hasDew = true;
                point.dewPressurePa = refineTransitionPressure(
                    z,
                    temperatureK,
                    pressures[i],
                    pressures[i + 1],
                    options.vaporPhaseCode,
                    options.twoPhaseCode,
                    options);
                ++point.boundedTwoPhaseIntervals;
            }
            if (detail::phaseMatches(codeA, options.twoPhaseCode) && detail::phaseMatches(codeB, options.liquidPhaseCode))
            {
                point.hasBubble = true;
                point.bubblePressurePa = refineTransitionPressure(
                    z,
                    temperatureK,
                    pressures[i],
                    pressures[i + 1],
                    options.liquidPhaseCode,
                    options.twoPhaseCode,
                    options);
            }
        }
        return point;
    }

    [[nodiscard]] double refinePresenceBitTransition(
        const Composition &z,
        double temperatureK,
        double lowerPressure,
        double upperPressure,
        std::uint8_t phaseBit,
        const MultiphaseBoundaryOptions &options) const
    {
        double lower = lowerPressure;
        double upper = upperPressure;
        Sample lowerSample = evaluate(lower, temperatureK, z);
        Sample upperSample = evaluate(upper, temperatureK, z);
        if (lowerSample.statusCode != 0 || upperSample.statusCode != 0)
            return 0.5 * (lower + upper);

        const auto present = [phaseBit](const Sample &sample) {
            return (static_cast<std::uint8_t>(sample.phaseCode()) & phaseBit) != 0u;
        };
        const bool lowerPresence = present(lowerSample);
        const bool upperPresence = present(upperSample);
        if (lowerPresence == upperPresence)
            return 0.5 * (lower + upper);

        for (int iter = 0; iter < options.maxRefinementIterations; ++iter)
        {
            if (detail::intervalSmallEnough(lower, upper, options.relativePressureTolerance))
                break;
            const double middle = 0.5 * (lower + upper);
            const Sample middleSample = evaluate(middle, temperatureK, z);
            if (middleSample.statusCode != 0)
                break;
            if (present(middleSample) == lowerPresence)
                lower = middle;
            else
                upper = middle;
        }
        return 0.5 * (lower + upper);
    }

    [[nodiscard]] double refineTransitionPressure(
        const Composition &z,
        double temperatureK,
        double lowerPressure,
        double upperPressure,
        int singlePhaseCode,
        int twoPhaseCode,
        const EnvelopeOptions &options) const
    {
        double pSingle = lowerPressure;
        double pTwo = upperPressure;
        Sample lower = evaluate(lowerPressure, temperatureK, z);
        Sample upper = evaluate(upperPressure, temperatureK, z);

        if (!(lower.statusCode == 0 && upper.statusCode == 0))
            return 0.5 * (lowerPressure + upperPressure);

        if (lower.phaseCode() == twoPhaseCode && upper.phaseCode() == singlePhaseCode)
        {
            pSingle = upperPressure;
            pTwo = lowerPressure;
            std::swap(lower, upper);
        }

        if (!(lower.phaseCode() == singlePhaseCode && upper.phaseCode() == twoPhaseCode))
            return 0.5 * (lowerPressure + upperPressure);

        for (int iter = 0; iter < options.maxRefinementIterations; ++iter)
        {
            if (detail::intervalSmallEnough(pSingle, pTwo, options.relativePressureTolerance))
                break;
            const double middlePressure = 0.5 * (pSingle + pTwo);
            const Sample middle = evaluate(middlePressure, temperatureK, z);
            if (middle.statusCode != 0)
                break;

            if (middle.phaseCode() == twoPhaseCode)
                pTwo = middlePressure;
            else if (middle.phaseCode() == singlePhaseCode)
                pSingle = middlePressure;
            else
                break;
        }
        return 0.5 * (pSingle + pTwo);
    }

    const Flash &flash_;
    PhaseDiagramFlashPolicy policy_{};
};

/** @brief 规则 P-T 相图的 CSV 写出器。 */
template <class Indices, class Name>
inline void writeCsv(
    const std::filesystem::path &path,
    const PressureTemperatureMap<Indices> &map,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create phase-diagram CSV: " + path.string());
    out << std::setprecision(17);
    out << "temperature_index,pressure_index,";
    detail::writeCommonHeader<Indices>(out, componentNames);
    out << '\n';

    const std::size_t nPressure = map.pressureAxis.points;
    for (std::size_t row = 0; row < map.samples.size(); ++row)
    {
        out << row / nPressure << ',' << row % nPressure << ',';
        detail::writeCommonRow(out, map.samples[row]);
        out << '\n';
    }
}

/** @brief 压力-组成扫描的 CSV 写出器。 */
template <class Indices, class Name>
inline void writeCsv(
    const std::filesystem::path &path,
    const PressureCompositionMap<Indices> &map,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create phase-diagram CSV: " + path.string());
    out << std::setprecision(17);
    out << "path_index,pressure_index,path_fraction,";
    detail::writeCommonHeader<Indices>(out, componentNames);
    out << '\n';

    const std::size_t nPressure = map.pressureAxis.points;
    for (std::size_t row = 0; row < map.samples.size(); ++row)
    {
        out << row / nPressure << ',' << row % nPressure << ',' << map.pathFraction[row] << ',';
        detail::writeCommonRow(out, map.samples[row]);
        out << '\n';
    }
}

/** @brief 固定 P、T 三元组成相图的 CSV 写出器。 */
template <class Indices, class Name>
inline void writeCsv(
    const std::filesystem::path &path,
    const TernaryCompositionMap<Indices> &map,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create phase-diagram CSV: " + path.string());
    out << std::setprecision(17);
    out << "point_index,ternary_a,ternary_b,ternary_c,ternary_x,ternary_y,";
    detail::writeCommonHeader<Indices>(out, componentNames);
    out << '\n';

    constexpr double rootThreeOverTwo = 0.86602540378443864676;
    for (std::size_t row = 0; row < map.samples.size(); ++row)
    {
        const auto fraction = map.ternaryFraction[row];
        const double x = fraction[1] + 0.5 * fraction[2];
        const double y = rootThreeOverTwo * fraction[2];
        out << row << ',' << fraction[0] << ',' << fraction[1] << ',' << fraction[2]
            << ',' << x << ',' << y << ',';
        detail::writeCommonRow(out, map.samples[row]);
        out << '\n';
    }
}

/** @brief 固定总体组成泡点/露点曲线的 CSV 写出器。 */
template <class Indices, class Name>
inline void writeCsv(
    const std::filesystem::path &path,
    const PressureTemperatureEnvelope<Indices> &envelope,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create envelope CSV: " + path.string());
    out << std::setprecision(17);
    out << "temperature_index,temperature_K,has_dew,dew_pressure_Pa,dew_pressure_bar"
        << ",has_bubble,bubble_pressure_Pa,bubble_pressure_bar"
        << ",envelope_width_Pa,envelope_width_bar,bounded_two_phase_intervals"
        << ",critical_point,critical_temperature_K,critical_pressure_Pa,critical_pressure_bar"
        << ",critical_envelope_width_Pa,critical_envelope_width_bar"
        << ",critical_relative_envelope_width,critical_oil_gas_composition_L1,critical_estimate_score,critical_estimate_within_tolerance";
    for (const auto &name : componentNames)
        out << ",z_" << detail::sanitizeColumnName(name);
    out << '\n';

    for (std::size_t i = 0; i < envelope.points.size(); ++i)
    {
        const auto &p = envelope.points[i];
        const bool isCritical = envelope.hasCriticalPoint &&
            std::abs(p.coordinate - envelope.criticalTemperatureK) <= 1.0e-12 * std::max(1.0, std::abs(envelope.criticalTemperatureK));
        out << i << ',' << p.coordinate << ',' << (p.hasDew ? 1 : 0) << ',' << p.dewPressurePa << ',' << p.dewPressurePa / 1.0e5
            << ',' << (p.hasBubble ? 1 : 0) << ',' << p.bubblePressurePa << ',' << p.bubblePressurePa / 1.0e5
            << ',' << p.envelopeWidthPa() << ',' << p.envelopeWidthPa() / 1.0e5
            << ',' << p.boundedTwoPhaseIntervals
            << ',' << (isCritical ? 1 : 0) << ',' << envelope.criticalTemperatureK
            << ',' << envelope.criticalPressurePa << ',' << envelope.criticalPressurePa / 1.0e5
            << ',' << envelope.criticalEnvelopeWidthPa << ',' << envelope.criticalEnvelopeWidthPa / 1.0e5
            << ',' << envelope.criticalRelativeEnvelopeWidth
            << ',' << envelope.criticalOilGasCompositionL1
            << ',' << envelope.criticalEstimateScore
            << ',' << (envelope.criticalEstimateWithinTolerance ? 1 : 0);
        for (double value : envelope.overallComposition)
            out << ',' << value;
        out << '\n';
    }
}

/** @brief 组成路径泡点/露点曲线的 CSV 写出器。 */
template <class Indices, class Name>
inline void writeCsv(
    const std::filesystem::path &path,
    const PressureCompositionEnvelope<Indices> &envelope,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create composition-envelope CSV: " + path.string());
    out << std::setprecision(17);
    out << "path_index,path_fraction,temperature_K,has_dew,dew_pressure_Pa,dew_pressure_bar"
        << ",has_bubble,bubble_pressure_Pa,bubble_pressure_bar,envelope_width_Pa,envelope_width_bar"
        << ",bounded_two_phase_intervals";
    for (const auto &name : componentNames)
        out << ",z_" << detail::sanitizeColumnName(name);
    out << '\n';

    const auto fractions = envelope.pathAxis.values();
    for (std::size_t i = 0; i < envelope.points.size(); ++i)
    {
        const auto &p = envelope.points[i];
        const auto z = detail::interpolatedComposition(envelope.startComposition, envelope.endComposition, fractions[i]);
        out << i << ',' << fractions[i] << ',' << envelope.temperatureK
            << ',' << (p.hasDew ? 1 : 0) << ',' << p.dewPressurePa << ',' << p.dewPressurePa / 1.0e5
            << ',' << (p.hasBubble ? 1 : 0) << ',' << p.bubblePressurePa << ',' << p.bubblePressurePa / 1.0e5
            << ',' << p.envelopeWidthPa() << ',' << p.envelopeWidthPa() / 1.0e5
            << ',' << p.boundedTwoPhaseIntervals;
        for (double value : z)
            out << ',' << value;
        out << '\n';
    }
}

/** @brief 组成路径近似临界轨迹的 CSV 写出器。 */
template <class Indices, class Name>
inline void writeCsv(
    const std::filesystem::path &path,
    const CriticalLocusPath<Indices> &locus,
    const std::array<Name, Indices::numComponents> &componentNames)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create critical-locus CSV: " + path.string());
    out << std::setprecision(17);
    out << "path_index,path_fraction,valid,critical_temperature_K,critical_pressure_Pa,critical_pressure_bar"
        << ",critical_envelope_width_Pa,critical_envelope_width_bar"
        << ",critical_relative_envelope_width,critical_oil_gas_composition_L1,critical_estimate_score,critical_estimate_within_tolerance";
    for (const auto &name : componentNames)
        out << ",z_" << detail::sanitizeColumnName(name);
    out << '\n';

    for (std::size_t i = 0; i < locus.points.size(); ++i)
    {
        const auto &p = locus.points[i];
        out << i << ',' << p.pathFraction << ',' << (p.valid ? 1 : 0)
            << ',' << p.criticalTemperatureK << ',' << p.criticalPressurePa << ',' << p.criticalPressurePa / 1.0e5
            << ',' << p.criticalEnvelopeWidthPa << ',' << p.criticalEnvelopeWidthPa / 1.0e5
            << ',' << p.criticalRelativeEnvelopeWidth << ',' << p.criticalOilGasCompositionL1
            << ',' << p.criticalEstimateScore
            << ',' << (p.criticalEstimateWithinTolerance ? 1 : 0);
        for (double value : p.overallComposition)
            out << ',' << value;
        out << '\n';
    }
}


/** @brief 非限制 O/G/W 相出现边界的 CSV 写出器。 */
template <class Indices>
inline void writeCsv(
    const std::filesystem::path &path,
    const PressureTemperaturePhaseBoundaries<Indices> &boundaries)
{
    detail::prepareOutputFile(path);
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create multiphase-boundary CSV: " + path.string());
    out << std::setprecision(17);
    out << "temperature_K,pressure_Pa,pressure_bar,boundary_code,boundary_name,lower_phase_code,upper_phase_code\n";
    for (const auto &point : boundaries.points)
    {
        const char *name = point.kind == MultiphaseBoundaryKind::OilOnset ? "oil_onset" :
            (point.kind == MultiphaseBoundaryKind::GasOnset ? "gas_onset" : "water_onset");
        out << point.temperatureK << ',' << point.pressurePa << ',' << point.pressurePa / 1.0e5
            << ',' << static_cast<int>(point.kind) << ',' << name << ','
            << point.lowerPhaseCode << ',' << point.upperPhaseCode << '\n';
    }
}

} // namespace MPMC::tools
