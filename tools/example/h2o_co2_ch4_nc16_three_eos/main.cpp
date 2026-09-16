/**
 * @file main.cpp
 * @brief H2O–CO2–CH4–nC16 的 PR/SW/CPA 三 EOS 对比扫描入口。
 */
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <tools/matlab_phase_diagram.hpp>
#include <tools/phase_diagram.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>

namespace
{

using Config = MPMC::CompositionalModelConfig<
    4,
    true,
    false,
    false,
    false,
    false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, Indices::numComponents>;

constexpr std::size_t N = static_cast<std::size_t>(Indices::numComponents);
inline const std::array<std::string, N> componentNames{"H2O", "CO2", "CH4", "nC16"};

inline constexpr std::array<double, N> tc{647.096, 304.1282, 190.564, 723.0};
inline constexpr std::array<double, N> pc{22.064e6, 7.3773e6, 4.5992e6, 1.410e6};
inline constexpr std::array<double, N> vc{5.6e-5, 9.4e-5, 9.9e-5, 9.0e-4};
inline constexpr std::array<double, N> omega{0.3443, 0.22394, 0.01142, 0.742};
inline constexpr std::array<double, N> mw{0.01801528, 0.0440095, 0.016043, 0.226441};
inline constexpr std::array<std::array<double, N>, N> kij{{
    {{0.0, 0.1896, 0.4850, 0.5000}},
    {{0.1896, 0.0, 0.1200, 0.0900}},
    {{0.4850, 0.1200, 0.0, 0.0000}},
    {{0.5000, 0.0900, 0.0000, 0.0}}
}};

MPMC::CompositionalMixture<Indices> makeMixture()
{
    return MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij);
}

Eos makePrEos()
{
    return Eos(
        0.45724,
        0.07780,
        makeMixture(),
        5,
        2.414213562373095,
        -0.414213562373095,
        1.0e-30);
}

Eos makeSwEos()
{
    Eos eos = makePrEos();
    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    sw.aqueousWaterBip[1] = [](double temperature, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
            temperature, tc[1], salinity);
    };
    sw.aqueousWaterBip[2] = [](double temperature, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperature, tc[2], omega[2], salinity);
    };
    // The published SW hydrocarbon correlation was fitted through n-butane.
    // Avoid silently extrapolating it to nC16: use the same transparent heavy
    // water/hydrocarbon engineering BIP as the ordinary-PR benchmark instead.
    sw.aqueousWaterBip[3] = [](double, double) { return 0.5000; };
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

Eos makeCpaEos()
{
    Eos eos = makePrEos();
    Eos::CubicPlusAssociationOptions cpa;
    constexpr double R = MPMC::units::gasConstant;
    for (std::size_t c = 0; c < N; ++c)
    {
        cpa.a0[c] = 0.42748 * R * R * tc[c] * tc[c] / pc[c];
        cpa.b[c] = 0.08664 * R * tc[c] / pc[c];
        cpa.c1[c] = 0.480 + 1.574 * omega[c] - 0.176 * omega[c] * omega[c];
    }
    cpa.a0[0] = MPMC::StandardCpaWater4C::a0;
    cpa.b[0] = MPMC::StandardCpaWater4C::b;
    cpa.c1[0] = MPMC::StandardCpaWater4C::c1;
    cpa.associationEnergy[0] = MPMC::StandardCpaWater4C::epsilon;
    cpa.associationVolume[0] = MPMC::StandardCpaWater4C::beta;
    cpa.donorSites[0] = MPMC::StandardCpaWater4C::donorSites;
    cpa.acceptorSites[0] = MPMC::StandardCpaWater4C::acceptorSites;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

std::size_t countFailures(const MPMC::tools::PressureTemperatureMap<Indices> &map)
{
    std::size_t result = 0;
    for (const auto &sample : map.samples)
        result += sample.statusCode == 0 ? 0u : 1u;
    return result;
}

std::size_t countPhase(const MPMC::tools::PressureTemperatureMap<Indices> &map, int code)
{
    std::size_t result = 0;
    for (const auto &sample : map.samples)
        result += sample.phaseCode() == code ? 1u : 0u;
    return result;
}

struct BackendSummary
{
    std::string name;
    std::string directory;
    std::size_t ptSamples{0};
    std::size_t failures{0};
    std::size_t threePhase{0};
    std::size_t waterBoundaries{0};
    bool criticalEstimate{false};
};

template <class EosFactory>
BackendSummary runBackend(
    std::string_view shortName,
    std::string_view displayName,
    EosFactory &&factory,
    const std::filesystem::path &root)
{
    const std::filesystem::path outDir = root / std::string(shortName);
    std::filesystem::create_directories(outDir);

    const Eos eos = factory();
    MPMC::ThreePhaseFlashOptions flashOptions;
    flashOptions.waterComponent = 0;
    const Flash flash(eos, flashOptions);

    MPMC::tools::PhaseDiagramSampler<Indices> threePhaseSampler(
        flash,
        MPMC::tools::PhaseDiagramFlashPolicy::unrestricted());
    const auto oilGas = MPMC::PhasePresence(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit);
    MPMC::tools::PhaseDiagramSampler<Indices> vleSampler(
        flash,
        MPMC::tools::PhaseDiagramFlashPolicy::restricted(oilGas));
    MPMC::tools::EnvelopeOptions envelopeOptions;
    envelopeOptions.maxRefinementIterations = 20;
    envelopeOptions.relativePressureTolerance = 1.0e-5;

    // Representative wet reservoir fluid. The grid deliberately extends above
    // the ordinary-water critical T/P reference in both coordinates.
    const Composition z{0.10, 0.25, 0.50, 0.15};
    const MPMC::tools::ScanAxis broadT{300.0, 850.0, 11, MPMC::tools::AxisSpacing::Linear};
    const MPMC::tools::ScanAxis broadP{1.0e5, 6.0e7, 13, MPMC::tools::AxisSpacing::Logarithmic};

    const auto pt = threePhaseSampler.pressureTemperature(z, broadT, broadP);
    MPMC::tools::writeCsv(outDir / "pt_three_phase.csv", pt, componentNames);
    MPMC::tools::writePressureTemperatureMatlab(
        outDir / "plot_pt_three_phase.m",
        "pt_three_phase.csv",
        std::string(displayName) + " unrestricted O/G/W P-T map");
    MPMC::tools::writePressureTemperatureQualityMatlab(
        outDir / "plot_pt_phase_count_quality.m",
        "pt_three_phase.csv",
        std::string(displayName) + " phase count and gas-quality contours",
        true);

    MPMC::tools::MultiphaseBoundaryOptions boundaryOptions;
    boundaryOptions.maxRefinementIterations = 1;
    boundaryOptions.relativePressureTolerance = 1.0e-3;
    const auto boundaries = threePhaseSampler.pressureTemperaturePhaseBoundaries(pt, boundaryOptions);
    MPMC::tools::writeCsv(outDir / "pt_phase_onset_boundaries.csv", boundaries);
    MPMC::tools::writeMultiphaseBoundariesMatlab(
        outDir / "plot_pt_phase_onset_boundaries.m",
        "pt_phase_onset_boundaries.csv",
        std::string(displayName) + " O/G/W onset boundaries");

    // Standard paper-style bubble/dew envelope: the same four-component EOS,
    // but projected onto O+G VLE so a water-rich third phase does not redefine
    // the conventional bubble/dew terminology.
    const MPMC::tools::ScanAxis envelopeT{280.0, 800.0, 15, MPMC::tools::AxisSpacing::Linear};
    const MPMC::tools::ScanAxis envelopeP{1.0e5, 6.0e7, 17, MPMC::tools::AxisSpacing::Logarithmic};
    const auto vleMap = vleSampler.pressureTemperature(z, envelopeT, envelopeP);
    const auto envelope = vleSampler.pressureTemperatureEnvelope(vleMap, envelopeOptions);
    MPMC::tools::writeCsv(outDir / "pt_vle_envelope.csv", envelope, componentNames);
    MPMC::tools::writePressureTemperatureEnvelopeMatlab(
        outDir / "plot_pt_vle_envelope.m",
        "pt_vle_envelope.csv",
        std::string(displayName) + " O/G bubble-dew envelope projection",
        true);

    // Two common flash-result curves: pressure sweep above water Tc and a
    // temperature sweep at 300 bar.  CSV includes beta, S, rho_m, Z, K and all
    // phase compositions for direct paper plotting.
    const auto pressureProfile = threePhaseSampler.pressureTemperature(
        z,
        MPMC::tools::ScanAxis{673.15, 673.15, 1, MPMC::tools::AxisSpacing::Linear},
        MPMC::tools::ScanAxis{1.0e5, 6.0e7, 21, MPMC::tools::AxisSpacing::Logarithmic});
    MPMC::tools::writeCsv(outDir / "pressure_profile_673K.csv", pressureProfile, componentNames);
    MPMC::tools::writePressureFlashProfileMatlab(
        outDir / "plot_pressure_profile_673K.m",
        "pressure_profile_673K.csv",
        componentNames,
        std::string(displayName) + " O/G/W flash at 673.15 K");

    const auto temperatureProfile = threePhaseSampler.pressureTemperature(
        z,
        MPMC::tools::ScanAxis{300.0, 850.0, 25, MPMC::tools::AxisSpacing::Linear},
        MPMC::tools::ScanAxis{3.0e7, 3.0e7, 1, MPMC::tools::AxisSpacing::Linear});
    MPMC::tools::writeCsv(outDir / "temperature_profile_300bar.csv", temperatureProfile, componentNames);
    MPMC::tools::writeTemperatureFlashProfileMatlab(
        outDir / "plot_temperature_profile_300bar.m",
        "temperature_profile_300bar.csv",
        componentNames,
        std::string(displayName) + " O/G/W flash at 300 bar");

    // Approximate critical-locus trend along a wet-gas -> wet-oil path.  This
    // remains explicitly an envelope-coalescence estimate, not a rigorous
    // Michelsen criticality-equation solver.
    const Composition wetGas{0.05, 0.25, 0.60, 0.10};
    const Composition wetOil{0.05, 0.20, 0.20, 0.55};
    const auto critical = vleSampler.criticalLocus(
        MPMC::tools::ScanAxis{0.0, 1.0, 3, MPMC::tools::AxisSpacing::Linear},
        wetGas,
        wetOil,
        MPMC::tools::ScanAxis{280.0, 800.0, 5, MPMC::tools::AxisSpacing::Linear},
        MPMC::tools::ScanAxis{1.0e5, 6.0e7, 7, MPMC::tools::AxisSpacing::Logarithmic},
        envelopeOptions);
    MPMC::tools::writeCsv(outDir / "critical_locus_estimate.csv", critical, componentNames);
    MPMC::tools::writeCriticalLocusMatlab(
        outDir / "plot_critical_locus_estimate.m",
        "critical_locus_estimate.csv",
        std::string(displayName) + " critical-locus estimate",
        true);

    std::size_t waterBoundaries = 0;
    for (const auto &point : boundaries.points)
        waterBoundaries += point.kind == MPMC::tools::MultiphaseBoundaryKind::WaterOnset ? 1u : 0u;

    BackendSummary summary;
    summary.name = std::string(displayName);
    summary.directory = std::string(shortName);
    summary.ptSamples = pt.samples.size();
    summary.failures = countFailures(pt);
    summary.threePhase = countPhase(pt, 7);
    summary.waterBoundaries = waterBoundaries;
    summary.criticalEstimate = envelope.hasCriticalPoint;
    std::cout << "[three-eos-phase-tools] completed " << displayName
              << ": P-T samples=" << summary.ptSamples
              << ", failures=" << summary.failures
              << ", O+G+W=" << summary.threePhase
              << ", water-onset points=" << summary.waterBoundaries
              << ", critical-estimate=" << (summary.criticalEstimate ? 1 : 0) << std::endl;
    return summary;
}

} // namespace

int main(int argc, char **argv)
{
    const std::filesystem::path outputDirectory =
        argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("./results");
    std::filesystem::create_directories(outputDirectory);

    const std::array<BackendSummary, 3> summaries{
        runBackend("pr", "Peng-Robinson", makePrEos, outputDirectory),
        runBackend("sw", "Soreide-Whitson", makeSwEos, outputDirectory),
        runBackend("cpa", "Cubic-Plus-Association", makeCpaEos, outputDirectory)};

    for (const auto &s : summaries)
    {
        std::cout << "[three-eos-phase-tools] " << s.name
                  << ": P-T samples=" << s.ptSamples
                  << ", failures=" << s.failures
                  << ", O+G+W=" << s.threePhase
                  << ", water-onset points=" << s.waterBoundaries
                  << ", critical-estimate=" << (s.criticalEstimate ? 1 : 0) << '\n';
    }

    const std::array<std::string, 3> envelopeFiles{
        "pr/pt_vle_envelope.csv", "sw/pt_vle_envelope.csv", "cpa/pt_vle_envelope.csv"};
    const std::array<std::string, 3> criticalFiles{
        "pr/critical_locus_estimate.csv", "sw/critical_locus_estimate.csv", "cpa/critical_locus_estimate.csv"};
    const std::array<std::string, 3> backendNames{"PR", "SW", "CPA"};
    MPMC::tools::writeBackendEnvelopeComparisonMatlab(
        outputDirectory / "plot_compare_vle_envelopes.m",
        envelopeFiles,
        backendNames,
        "PR / Soreide-Whitson / CPA bubble-dew envelope comparison");
    MPMC::tools::writeCriticalLocusComparisonMatlab(
        outputDirectory / "plot_compare_critical_loci.m",
        criticalFiles,
        backendNames,
        "PR / Soreide-Whitson / CPA critical-locus estimates");

    std::cout << "[three-eos-phase-tools] H2O critical reference: T="
              << MPMC::tools::WaterCriticalReference::temperatureK
              << " K, P=" << MPMC::tools::WaterCriticalReference::pressureBar
              << " bar\n";
    std::cout << "[three-eos-phase-tools] generated under " << outputDirectory << '\n';
    return 0;
}
