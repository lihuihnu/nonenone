/**
 * @file main.cpp
 * @brief C2–C3–nC5 PR 相图与 flash 验证示例入口。
 */
#include <indices/indices.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <tools/matlab_phase_diagram.hpp>
#include <tools/phase_diagram.hpp>

#include <array>
#include <filesystem>
#include <iostream>
#include <string>

namespace
{

using Config = MPMC::CompositionalModelConfig<3, false, false>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, Indices::numComponents>;

inline const std::array<std::string, 3> componentNames{
    "C2", "C3", "nC5"};

Eos makeEos()
{
    // Simple Peng-Robinson hydrocarbon set also used by the public
    // Petroleum Office regression in the project test suite.
    const std::array<double, 3> tc{305.32, 369.83, 469.70};
    const std::array<double, 3> pc{4.872e6, 4.248e6, 3.370e6};
    const std::array<double, 3> vc{1.458e-4, 2.000e-4, 3.110e-4};
    const std::array<double, 3> omega{0.099, 0.152, 0.252};
    const std::array<double, 3> mw{0.03007, 0.04410, 0.07215};
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

void printSummary(const char *name, std::size_t samples, std::size_t failures)
{
    std::cout << "[phase-diagram] " << name << ": " << samples
              << " samples, failures=" << failures << '\n';
}

template <class Map>
std::size_t countFailures(const Map &map)
{
    std::size_t failures = 0;
    for (const auto &sample : map.samples)
        failures += sample.statusCode == 0 ? 0u : 1u;
    return failures;
}

} // namespace

int main(int argc, char **argv)
{
    const std::filesystem::path outputDirectory =
        argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("./results");
    std::filesystem::create_directories(outputDirectory);

    const Eos eos = makeEos();
    const Flash flash(eos);

    // This example is a conventional hydrocarbon vapor-liquid diagram.  Limit
    // the active phase set to O+G so the third liquid slot used by fully
    // compositional H2O systems cannot be interpreted as an extra phase.
    const auto oilGas = MPMC::PhasePresence(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit);
    MPMC::tools::PhaseDiagramSampler<Indices> sampler(
        flash,
        MPMC::tools::PhaseDiagramFlashPolicy::restricted(oilGas));
    const MPMC::tools::EnvelopeOptions envelopeOptions{};

    // 1. P-T map at fixed overall composition z=(0.3,0.4,0.3).
    const Composition z{0.30, 0.40, 0.30};
    const auto pt = sampler.pressureTemperature(
        z,
        MPMC::tools::ScanAxis{270.0, 330.0, 31, MPMC::tools::AxisSpacing::Linear},
        MPMC::tools::ScanAxis{1.0e5, 2.0e6, 41, MPMC::tools::AxisSpacing::Linear});
    MPMC::tools::writeCsv(outputDirectory / "pt_phase_map.csv", pt, componentNames);
    MPMC::tools::writePressureTemperatureMatlab(
        outputDirectory / "plot_pt_phase_map.m",
        "pt_phase_map.csv",
        "C2-C3-nC5 Peng-Robinson P-T phase map");
    printSummary("P-T", pt.samples.size(), countFailures(pt));

    // 1b. Continuous bubble/dew envelope and approximate critical point.
    const auto ptEnvelope = sampler.pressureTemperatureEnvelope(pt, envelopeOptions);
    MPMC::tools::writeCsv(outputDirectory / "pt_envelope.csv", ptEnvelope, componentNames);
    MPMC::tools::writePressureTemperatureEnvelopeMatlab(
        outputDirectory / "plot_pt_envelope.m",
        "pt_envelope.csv",
        "C2-C3-nC5 Peng-Robinson bubble/dew envelope");
    std::cout << "[phase-diagram] P-T envelope: " << ptEnvelope.points.size()
              << " temperatures, critical_detected=" << (ptEnvelope.hasCriticalPoint ? 1 : 0) << '\n';

    // 2. Pressure-composition map along the C2<->nC5 binary edge at 320 K.
    const Composition c2Rich{0.999, 0.0, 0.001};
    const Composition nc5Rich{0.001, 0.0, 0.999};
    const auto px = sampler.pressureComposition(
        320.0,
        MPMC::tools::ScanAxis{1.0e5, 5.0e6, 51, MPMC::tools::AxisSpacing::Linear},
        c2Rich,
        nc5Rich,
        MPMC::tools::ScanAxis{0.0, 1.0, 51, MPMC::tools::AxisSpacing::Linear});
    MPMC::tools::writeCsv(outputDirectory / "px_phase_map.csv", px, componentNames);
    MPMC::tools::writePressureCompositionMatlab(
        outputDirectory / "plot_px_phase_map.m",
        "px_phase_map.csv",
        "path fraction toward nC5-rich endpoint",
        "C2-nC5 Peng-Robinson P-composition map at 320 K");
    printSummary("P-composition", px.samples.size(), countFailures(px));

    // 2b. Continuous bubble/dew lines along the composition path.
    const auto pxEnvelope = sampler.pressureCompositionEnvelope(px, envelopeOptions);
    MPMC::tools::writeCsv(outputDirectory / "px_envelope.csv", pxEnvelope, componentNames);
    MPMC::tools::writePressureCompositionEnvelopeMatlab(
        outputDirectory / "plot_px_envelope.m",
        "px_envelope.csv",
        "path fraction toward nC5-rich endpoint",
        "C2-nC5 Peng-Robinson bubble/dew lines at 320 K");

    // 3. Ternary map at the same state as the public 8 bar, 320 K flash check.
    const auto ternary = sampler.ternaryComposition(
        8.0e5,
        320.0,
        {0, 1, 2},
        32);
    MPMC::tools::writeCsv(
        outputDirectory / "ternary_phase_map.csv", ternary, componentNames);
    MPMC::tools::writeTernaryMatlab(
        outputDirectory / "plot_ternary_phase_map.m",
        "ternary_phase_map.csv",
        componentNames,
        {0, 1, 2},
        "C2-C3-nC5 Peng-Robinson ternary phase map at 320 K and 8 bar");
    printSummary("ternary", ternary.samples.size(), countFailures(ternary));

    // 4. Approximate critical locus along the C2<->nC5 edge.
    const auto criticalLocus = sampler.criticalLocus(
        MPMC::tools::ScanAxis{0.0, 1.0, 11, MPMC::tools::AxisSpacing::Linear},
        c2Rich,
        nc5Rich,
        MPMC::tools::ScanAxis{270.0, 460.0, 31, MPMC::tools::AxisSpacing::Linear},
        MPMC::tools::ScanAxis{1.0e5, 5.0e6, 61, MPMC::tools::AxisSpacing::Linear},
        envelopeOptions);
    MPMC::tools::writeCsv(outputDirectory / "critical_locus.csv", criticalLocus, componentNames);
    MPMC::tools::writeCriticalLocusMatlab(
        outputDirectory / "plot_critical_locus.m",
        "critical_locus.csv",
        "C2-nC5 Peng-Robinson approximate critical locus");
    std::size_t validCritical = 0;
    for (const auto &point : criticalLocus.points)
        validCritical += point.valid ? 1u : 0u;
    std::cout << "[phase-diagram] critical locus: " << criticalLocus.points.size()
              << " path points, valid=" << validCritical << '\n';

    std::cout << "MATLAB: cd('" << outputDirectory.string()
              << "'); run('plot_pt_phase_map.m'); run('plot_pt_envelope.m'); "
                 "run('plot_px_phase_map.m'); run('plot_px_envelope.m'); "
                 "run('plot_ternary_phase_map.m'); run('plot_critical_locus.m');\n";
    return 0;
}
