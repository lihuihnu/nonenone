/**
 * @file tools_phase_diagram_test.cpp
 * @brief 单元测试：验证 `tools_phase_diagram` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <tools/matlab_phase_diagram.hpp>
#include <tools/phase_diagram.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{

using Config = MPMC::CompositionalModelConfig<3, false, false>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, 3>;

Eos makeEos()
{
    const Composition tc{305.32, 369.83, 469.70};
    const Composition pc{4.872e6, 4.248e6, 3.370e6};
    const Composition vc{1.458e-4, 2.000e-4, 3.110e-4};
    const Composition omega{0.099, 0.152, 0.252};
    const Composition mw{0.03007, 0.04410, 0.07215};
    const std::array<std::array<double, 3>, 3> kij{};
    return Eos(
        0.45724,
        0.07780,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        1,
        2.414213562373095,
        -0.414213562373095,
        1.0e-30);
}

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::size_t lineCount(const std::filesystem::path &path)
{
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open generated test file");
    std::size_t count = 0;
    std::string line;
    while (std::getline(in, line))
        ++count;
    return count;
}

std::string readAll(const std::filesystem::path &path)
{
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open generated MATLAB file");
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

} // namespace

int main()
{
    const Eos eos = makeEos();
    const Flash flash(eos);
    const auto oilGas = MPMC::PhasePresence(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit);
    MPMC::tools::PhaseDiagramSampler<Indices> sampler(
        flash,
        MPMC::tools::PhaseDiagramFlashPolicy::restricted(oilGas));

    const Composition z{0.3, 0.4, 0.3};
    const auto point = sampler.evaluate(8.0e5, 320.0, z);
    require(point.statusCode == 0 && point.flash.converged,
            "phase-diagram production flash sample must converge");
    require(point.phaseCode() == 3,
            "public C2/C3/nC5 state must be O+G in the restricted diagram");
    require(std::abs(point.flash.phaseMoleFraction[0] - 0.2904665759) < 2.0e-6,
            "tools sampler must preserve the public PR liquid-fraction regression");
    require(std::abs(point.flash.phaseMoleFraction[1] - 0.7095334241) < 2.0e-6,
            "tools sampler must preserve the public PR vapor-fraction regression");

    const auto pt = sampler.pressureTemperature(
        z,
        {300.0, 320.0, 2, MPMC::tools::AxisSpacing::Linear},
        {5.0e5, 1.0e6, 3, MPMC::tools::AxisSpacing::Linear});
    require(pt.samples.size() == 6, "P-T scan size regression");

    const auto px = sampler.pressureComposition(
        320.0,
        {5.0e5, 1.0e6, 2, MPMC::tools::AxisSpacing::Linear},
        Composition{0.999, 0.0, 0.001},
        Composition{0.001, 0.0, 0.999},
        {0.0, 1.0, 3, MPMC::tools::AxisSpacing::Linear});
    require(px.samples.size() == 6 && px.pathFraction.size() == 6,
            "P-composition scan size regression");

    const auto ternary = sampler.ternaryComposition(8.0e5, 320.0, {0, 1, 2}, 4);
    require(ternary.samples.size() == 15 && ternary.ternaryFraction.size() == 15,
            "ternary triangular-grid point-count regression");

    // Envelope regression on a small but meaningful grid.
    const auto ptEnv = sampler.pressureTemperatureEnvelope(
        z,
        {280.0, 330.0, 6, MPMC::tools::AxisSpacing::Linear},
        {1.0e5, 2.0e6, 25, MPMC::tools::AxisSpacing::Linear});
    require(ptEnv.points.size() == 6, "P-T envelope point-count regression");
    std::size_t ptTransitions = 0;
    for (const auto &p : ptEnv.points)
        ptTransitions += (p.hasDew || p.hasBubble) ? 1u : 0u;
    require(ptTransitions > 0, "P-T envelope must detect at least one bubble/dew temperature");

    const auto boundaryMap = sampler.pressureTemperature(
        z,
        {280.0, 330.0, 6, MPMC::tools::AxisSpacing::Linear},
        {1.0e5, 2.0e6, 25, MPMC::tools::AxisSpacing::Linear});
    const auto boundaries = sampler.pressureTemperaturePhaseBoundaries(boundaryMap);
    require(!boundaries.points.empty(), "phase-onset boundary extraction must detect a phase-mask transition");
    require(std::abs(MPMC::tools::WaterCriticalReference::temperatureK - 647.096) < 1.0e-12,
            "IAPWS water critical temperature reference regression");
    require(std::abs(MPMC::tools::WaterCriticalReference::pressureBar - 220.64) < 1.0e-12,
            "IAPWS water critical pressure reference regression");

    const auto pxEnv = sampler.pressureCompositionEnvelope(
        320.0,
        {1.0e5, 5.0e6, 25, MPMC::tools::AxisSpacing::Linear},
        Composition{0.999, 0.0, 0.001},
        Composition{0.001, 0.0, 0.999},
        {0.0, 1.0, 5, MPMC::tools::AxisSpacing::Linear});
    require(pxEnv.points.size() == 5, "P-composition envelope point-count regression");
    std::size_t pxTransitions = 0;
    for (const auto &p : pxEnv.points)
        pxTransitions += (p.hasDew || p.hasBubble) ? 1u : 0u;
    require(pxTransitions > 0, "P-composition envelope must detect at least one bubble/dew composition");

    const auto criticalLocus = sampler.criticalLocus(
        {0.0, 1.0, 3, MPMC::tools::AxisSpacing::Linear},
        Composition{0.999, 0.0, 0.001},
        Composition{0.001, 0.0, 0.999},
        {280.0, 430.0, 6, MPMC::tools::AxisSpacing::Linear},
        {1.0e5, 5.0e6, 25, MPMC::tools::AxisSpacing::Linear});
    require(criticalLocus.points.size() == 3, "critical-locus point-count regression");
    std::size_t validCritical = 0;
    for (const auto &p : criticalLocus.points)
        validCritical += p.valid ? 1u : 0u;
    require(validCritical > 0, "critical locus must detect at least one valid approximate critical point");

    const auto temp = std::filesystem::temp_directory_path() / "mpmc_tools_phase_diagram_test";
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);
    const std::array<std::string, 3> names{"C2", "C3", "nC5"};

    MPMC::tools::writeCsv(temp / "pt.csv", pt, names);
    MPMC::tools::writeCsv(temp / "px.csv", px, names);
    MPMC::tools::writeCsv(temp / "ternary.csv", ternary, names);
    MPMC::tools::writeCsv(temp / "pt_env.csv", ptEnv, names);
    MPMC::tools::writeCsv(temp / "px_env.csv", pxEnv, names);
    MPMC::tools::writeCsv(temp / "critical_locus.csv", criticalLocus, names);
    MPMC::tools::writeCsv(temp / "boundaries.csv", boundaries);
    require(lineCount(temp / "pt.csv") == 7, "P-T CSV row count regression");
    require(lineCount(temp / "px.csv") == 7, "P-composition CSV row count regression");
    require(lineCount(temp / "ternary.csv") == 16, "ternary CSV row count regression");
    require(lineCount(temp / "pt_env.csv") == 7, "P-T envelope CSV row count regression");
    require(lineCount(temp / "px_env.csv") == 6, "P-composition envelope CSV row count regression");
    require(lineCount(temp / "critical_locus.csv") == 4, "critical locus CSV row count regression");
    require(lineCount(temp / "boundaries.csv") == boundaries.points.size() + 1,
            "phase-onset boundary CSV row count regression");
    require(readAll(temp / "pt.csv").find("Kg_C2") != std::string::npos,
            "common CSV must export gas/oil K-values");
    require(readAll(temp / "pt_env.csv").find("critical_estimate_score") != std::string::npos,
            "envelope CSV must export critical-estimate residual metrics");
    require(readAll(temp / "pt_env.csv").find("critical_estimate_within_tolerance") != std::string::npos,
            "envelope CSV must export critical-estimate quality flag");

    MPMC::tools::writePressureTemperatureMatlab(temp / "plot_pt.m", "pt.csv", "test PT");
    MPMC::tools::writePressureCompositionMatlab(temp / "plot_px.m", "px.csv", "lambda", "test PX");
    MPMC::tools::writePressureTemperatureEnvelopeMatlab(temp / "plot_pt_env.m", "pt_env.csv", "test envelope");
    MPMC::tools::writePressureCompositionEnvelopeMatlab(temp / "plot_px_env.m", "px_env.csv", "lambda", "test path envelope");
    MPMC::tools::writeCriticalLocusMatlab(temp / "plot_critical.m", "critical_locus.csv", "test critical locus");
    MPMC::tools::writePressureTemperatureQualityMatlab(temp / "plot_quality.m", "pt.csv", "test quality");
    MPMC::tools::writePressureTemperatureQualityMatlab(temp / "plot_quality_water.m", "pt.csv", "test quality water", true);
    MPMC::tools::writePressureFlashProfileMatlab(temp / "plot_pressure_profile.m", "pt.csv", names, "test pressure profile");
    MPMC::tools::writeTemperatureFlashProfileMatlab(temp / "plot_temperature_profile.m", "pt.csv", names, "test temperature profile");
    MPMC::tools::writeMultiphaseBoundariesMatlab(temp / "plot_boundaries.m", "boundaries.csv", "test boundaries");
    MPMC::tools::writeBackendEnvelopeComparisonMatlab(
        temp / "plot_compare_envelope.m",
        std::array<std::string, 2>{"pt_env.csv", "pt_env.csv"},
        std::array<std::string, 2>{"A", "B"},
        "test comparison");
    MPMC::tools::writeCriticalLocusComparisonMatlab(
        temp / "plot_compare_critical.m",
        std::array<std::string, 2>{"critical_locus.csv", "critical_locus.csv"},
        std::array<std::string, 2>{"A", "B"},
        "test critical comparison");
    MPMC::tools::writeTernaryMatlab(temp / "plot_ternary.m", "ternary.csv", names, {0, 1, 2});
    require(readAll(temp / "plot_pt.m").find("readtable('pt.csv')") != std::string::npos,
            "P-T MATLAB script must reference generated CSV");
    require(readAll(temp / "plot_pt_env.m").find("dew line") != std::string::npos,
            "P-T envelope MATLAB script must label dew line");
    require(readAll(temp / "plot_critical.m").find("critical locus") != std::string::npos,
            "critical-locus MATLAB script must label the critical locus");
    require(readAll(temp / "plot_ternary.m").find("data.xo_C2") != std::string::npos,
            "ternary MATLAB script must consume exported phase compositions");
    require(readAll(temp / "plot_quality.m").find("betaG") != std::string::npos,
            "quality MATLAB script must consume gas phase fraction");
    require(readAll(temp / "plot_quality.m").find("H2O critical reference") == std::string::npos,
            "generic MATLAB plots must not inject a water critical marker by default");
    require(readAll(temp / "plot_quality_water.m").find("H2O critical reference") != std::string::npos,
            "water-aware MATLAB plots must support an explicit water critical marker");
    require(readAll(temp / "plot_pressure_profile.m").find("data.Kg_C2") != std::string::npos,
            "pressure-profile MATLAB script must consume exported K-values");
    require(readAll(temp / "plot_boundaries.m").find("water onset") != std::string::npos,
            "boundary MATLAB script must expose water-onset plotting");
    require(readAll(temp / "plot_compare_envelope.m").find("critical estimate") != std::string::npos,
            "backend comparison script must label approximate critical states honestly");

    std::filesystem::remove_all(temp);
    std::cout << "tools_phase_diagram_test: PASS\n";
    return 0;
}
