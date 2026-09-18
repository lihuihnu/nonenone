/**
 * @file main.cpp
 * @brief External CPA proxy benchmark for Athabasca-bitumen + water.
 *
 * Reproduces Jia & Okuno (2018) Case 1 with the production CPA backend.
 * The literature parameters are deliberately isolated from the generated
 * kerogen Heavy pseudo-component; this executable is an implementation
 * benchmark, not a parameter-transfer or calibration tool.
 */
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Config = MPMC::CompositionalModelConfig<
    5, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Mixture = MPMC::CompositionalMixture<Indices>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, 5>;
using Matrix = std::array<std::array<double, 5>, 5>;

constexpr std::size_t water = 0;
constexpr std::size_t pc1 = 1;
constexpr std::size_t pc2 = 2;
constexpr std::size_t pc3 = 3;
constexpr std::size_t asphaltene = 4;

struct ReferencePoint
{
    double temperatureK{};
    double pressureMPa{};
    double experimentalWaterInBitumen{};
    double publishedCpaWaterInBitumen{};
};

struct BoundaryPoint
{
    double temperatureK{};
    double pressureMPa{};
    double temperatureUncertaintyK{};
    double pressureUncertaintyMPa{};
};

struct OwStabilityState
{
    bool valid{false};
    bool gasUnstable{false};
    double gasTrialSum{std::numeric_limits<double>::quiet_NaN()};
    double materialClosure{std::numeric_limits<double>::quiet_NaN()};
};

struct BoundaryPrediction
{
    bool found{false};
    double pressureMPa{std::numeric_limits<double>::quiet_NaN()};
    double bracketLowMPa{std::numeric_limits<double>::quiet_NaN()};
    double bracketHighMPa{std::numeric_limits<double>::quiet_NaN()};
    int transitionCount{0};
    int evaluations{0};
};

std::vector<std::string> splitCsv(const std::string &line)
{
    std::vector<std::string> fields;
    std::string current;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
        const char ch = line[i];
        if (ch == '"')
        {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"')
            {
                current.push_back('"');
                ++i;
            }
            else
                quoted = !quoted;
        }
        else if (ch == ',' && !quoted)
        {
            fields.push_back(current);
            current.clear();
        }
        else
            current.push_back(ch);
    }
    if (!current.empty() && current.back() == '\r')
        current.pop_back();
    fields.push_back(current);
    return fields;
}

std::vector<ReferencePoint> readReference(const std::filesystem::path &path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open Athabasca reference CSV: " + path.string());
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("Athabasca reference CSV is empty.");
    const auto header = splitCsv(line);
    auto column = [&](const std::string &name) {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end())
            throw std::runtime_error("Missing CSV column: " + name);
        return static_cast<std::size_t>(std::distance(header.begin(), it));
    };
    const auto cT = column("T_K");
    const auto cP = column("P_MPa");
    const auto cExp = column("xH2O_bitumen_rich_experiment");
    const auto cPub = column("xH2O_Jia2018_CPA");

    std::vector<ReferencePoint> points;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        const auto row = splitCsv(line);
        ReferencePoint p;
        p.temperatureK = std::stod(row.at(cT));
        p.pressureMPa = std::stod(row.at(cP));
        p.experimentalWaterInBitumen = std::stod(row.at(cExp));
        p.publishedCpaWaterInBitumen = std::stod(row.at(cPub));
        points.push_back(p);
    }
    return points;
}

std::vector<BoundaryPoint> readBoundaryReference(
    const std::filesystem::path &path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error(
            "Cannot open Athabasca Figure-7 boundary CSV: " + path.string());
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("Athabasca Figure-7 boundary CSV is empty.");
    const auto header = splitCsv(line);
    auto column = [&](const std::string &name) {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end())
            throw std::runtime_error("Missing boundary CSV column: " + name);
        return static_cast<std::size_t>(std::distance(header.begin(), it));
    };
    const auto cT = column("T_K");
    const auto cP = column("P_MPa");
    const auto cTu = column("T_uncertainty_K");
    const auto cPu = column("P_uncertainty_MPa");

    std::vector<BoundaryPoint> points;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        const auto row = splitCsv(line);
        BoundaryPoint p;
        p.temperatureK = std::stod(row.at(cT));
        p.pressureMPa = std::stod(row.at(cP));
        p.temperatureUncertaintyK = std::stod(row.at(cTu));
        p.pressureUncertaintyMPa = std::stod(row.at(cPu));
        points.push_back(p);
    }
    return points;
}

Composition bitumenComposition()
{
    Composition z{0.0, 0.2664, 0.4925, 0.1360, 0.1041};
    double sum = 0.0;
    for (std::size_t i = 1; i < z.size(); ++i)
        sum += z[i];
    for (std::size_t i = 1; i < z.size(); ++i)
        z[i] /= sum;
    return z;
}

Composition feedFromWaterMassFraction(double waterMassFraction)
{
    if (!(waterMassFraction > 0.0 && waterMassFraction < 1.0))
        throw std::invalid_argument(
            "Athabasca feed water mass fraction must lie in (0,1).");

    const Composition oil = bitumenComposition();
    constexpr std::array<double, 5> mwKgMol{
        0.01801528, 0.35243, 0.53903, 0.70684, 0.91619};
    double oilMw = 0.0;
    for (std::size_t i = 1; i < oil.size(); ++i)
        oilMw += oil[i] * mwKgMol[i];

    const double nWater = waterMassFraction / mwKgMol[0];
    const double nOil = (1.0 - waterMassFraction) / oilMw;
    const double xWater = nWater / (nWater + nOil);

    Composition z{};
    z[0] = xWater;
    for (std::size_t i = 1; i < z.size(); ++i)
        z[i] = (1.0 - xWater) * oil[i];
    return z;
}

Composition experimentalFeed()
{
    // Jia & Okuno Figure 7 caption states 55.9 wt% water + 44.1 wt%
    // Athabasca bitumen.  The Amani thesis Table 2.5 labels the same six
    // boundary points as 55.9 wt% bitumen + 44.1 wt% water.  Table-5 branch
    // parity retains the Jia convention; Figure-7 boundary parity audits both.
    return feedFromWaterMassFraction(0.559);
}

Eos makeJiaCase1Cpa()
{
    constexpr double R = MPMC::units::gasConstant;
    const std::array<double, 5> tc{
        647.096, 612.8, 745.2, 858.9, 1340.0};
    const std::array<double, 5> pc{
        22.064e6, 2.2871e6, 1.5885e6, 1.3327e6, 0.94189e6};
    // Vc is not used by the production CPA thermodynamic/flash equations.
    // Positive SRK-consistent critical-volume estimates are supplied only to
    // satisfy the generic mixture metadata contract.
    std::array<double, 5> vc{};
    for (std::size_t i = 0; i < vc.size(); ++i)
        vc[i] = R * tc[i] / (3.0 * pc[i]);
    const std::array<double, 5> omega{
        0.3443, 0.9673, 1.1751, 1.2018, 1.5567};
    const std::array<double, 5> mw{
        0.01801528, 0.35243, 0.53903, 0.70684, 0.91619};

    Matrix kij{};
    kij[water][pc1] = kij[pc1][water] = -0.0293;
    kij[water][pc2] = kij[pc2][water] = -0.0378;
    kij[water][pc3] = kij[pc3][water] = -0.0380;
    kij[water][asphaltene] = kij[asphaltene][water] = -0.0380;

    Eos eos(
        0.42748, 0.08664,
        Mixture(tc, pc, vc, omega, mw, kij),
        1, 1.0, 0.0, 1.0e-30);

    Eos::CubicPlusAssociationOptions cpa;
    for (std::size_t i = 0; i < 5; ++i)
    {
        cpa.a0[i] = 0.42748 * R * R * tc[i] * tc[i] / pc[i];
        cpa.b[i] = 0.08664 * R * tc[i] / pc[i];
        cpa.c1[i] =
            0.480 + 1.574 * omega[i] - 0.176 * omega[i] * omega[i];
    }

    // Jia & Okuno (2018) Case-1 water parameters as printed in the
    // article.  Keep these rounded publication values isolated in this
    // external reproduction benchmark rather than silently substituting the
    // repository's higher-precision StandardCpaWater4C constants.
    constexpr double jiaWaterA0 = 0.1227;   // 0.0001227 kPa m6/mol2
    constexpr double jiaWaterB = 1.45e-5;  // m3/mol
    constexpr double jiaWaterC1 = 0.6735;
    constexpr double jiaWaterEpsilon = 16655.0; // J/mol
    constexpr double jiaWaterBeta = 0.0692;
    cpa.a0[water] = jiaWaterA0;
    cpa.b[water] = jiaWaterB;
    cpa.c1[water] = jiaWaterC1;
    cpa.associationEnergy[water] = jiaWaterEpsilon;
    cpa.associationVolume[water] = jiaWaterBeta;
    cpa.donorSites[water] = 2;
    cpa.acceptorSites[water] = 2;

    // Jia & Okuno Table 4 asphaltene: 4C self association.
    cpa.a0[asphaltene] = 0.05202 * 1000.0; // kPa m6/mol2 -> Pa m6/mol2
    cpa.b[asphaltene] = 0.000914;
    cpa.c1[asphaltene] = 2.50;
    cpa.associationEnergy[asphaltene] = 26.0 * 1000.0; // kPa m3/mol -> J/mol
    cpa.associationVolume[asphaltene] = 0.05;
    cpa.donorSites[asphaltene] = 2;
    cpa.acceptorSites[asphaltene] = 2;

    // PC1-PC3 are inert (epsilon_AiBi=0) but solvating.  They carry one
    // electron-acceptor site so water donors can form the explicit cross bond.
    // Jia Table 4 gives beta_AiBi=0.07.  CR-1 then gives:
    // epsilon_water-PC = (epsilon_water + 0)/2;
    // beta_water-PC = sqrt(beta_water * 0.07).
    for (std::size_t pc : {pc1, pc2, pc3})
    {
        cpa.associationEnergy[pc] = 0.0;
        cpa.associationVolume[pc] = 0.07;
        cpa.donorSites[pc] = 0;
        cpa.acceptorSites[pc] = 1;
        cpa.crossAssociationEnergy[water][pc] =
            0.5 * jiaWaterEpsilon;
        cpa.crossAssociationVolume[water][pc] =
            std::sqrt(jiaWaterBeta * 0.07);
    }

    // Water/asphaltene cross association is intentionally left to the
    // production CR-1 fallback because both species are 4C self-associating.
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::SoaveRedlichKwong;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

double maxMaterialClosure(const Composition &z, const Flash::Result &result)
{
    double maximum = 0.0;
    for (std::size_t component = 0; component < z.size(); ++component)
    {
        double reconstructed = 0.0;
        for (std::size_t phase = 0; phase < 3; ++phase)
            reconstructed += result.phaseMoleFraction[phase] *
                result.composition[phase][component];
        maximum = std::max(
            maximum, std::abs(reconstructed - z[component]));
    }
    return maximum;
}

OwStabilityState evaluateOwStability(
    const Flash &flash,
    double pressureMPa,
    double temperatureK,
    const Composition &z)
{
    const auto oilWater = MPMC::PhasePresence(
        static_cast<std::uint8_t>(
            MPMC::PhasePresence::oilBit |
            MPMC::PhasePresence::waterBit));
    OwStabilityState state;
    const auto restricted = flash.flashRestricted(
        pressureMPa * 1.0e6, temperatureK, z, oilWater);
    if (!restricted.converged ||
        restricted.presence.bits() != oilWater.bits())
        return state;

    state.materialClosure = maxMaterialClosure(z, restricted);
    if (!std::isfinite(state.materialClosure) ||
        state.materialClosure > 1.0e-8)
        return state;

    const auto stability = flash.stabilityTest(
        pressureMPa * 1.0e6, temperatureK, z,
        restricted.presence, restricted.composition);
    if (!stability.valid)
        return state;

    const std::size_t gasSlot = static_cast<std::size_t>(
        MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));
    state.valid = true;
    state.gasUnstable = stability.missingPhaseUnstable[gasSlot];
    state.gasTrialSum = stability.trialSum[gasSlot];
    return state;
}

BoundaryPrediction findWlvWlBoundary(
    const Flash &flash,
    double temperatureK,
    const Composition &z)
{
    constexpr double pMinMPa = 2.0;
    constexpr double pMaxMPa = 30.0;
    constexpr double scanStepMPa = 0.10;
    constexpr int refinementLevels = 6;
    constexpr int refinementSegments = 20;

    BoundaryPrediction result;
    bool havePrevious = false;
    double previousPressure = 0.0;
    OwStabilityState previous;
    double selectedLow = std::numeric_limits<double>::quiet_NaN();
    double selectedHigh = std::numeric_limits<double>::quiet_NaN();

    // The experimental WLV-WL line is the high-pressure exit from the WLV
    // region.  CPA may contain additional lower-pressure stability flips, so
    // scan the full range and deliberately retain the HIGHEST-PRESSURE
    // gas-unstable -> gas-stable crossing rather than the first crossing.
    const int scanCount = static_cast<int>(
        std::llround((pMaxMPa - pMinMPa) / scanStepMPa));
    for (int k = 0; k <= scanCount; ++k)
    {
        const double pressure = pMinMPa + scanStepMPa * k;
        const auto state =
            evaluateOwStability(flash, pressure, temperatureK, z);
        ++result.evaluations;
        if (!state.valid)
            continue;

        if (havePrevious &&
            previous.gasUnstable && !state.gasUnstable)
        {
            ++result.transitionCount;
            selectedLow = previousPressure;
            selectedHigh = pressure;
        }
        previousPressure = pressure;
        previous = state;
        havePrevious = true;
    }

    if (result.transitionCount < 1 ||
        !std::isfinite(selectedLow) || !std::isfinite(selectedHigh))
        return result;

    // Refine by repeated local scans instead of pure bisection.  Near an
    // incipient phase the nonlinear restricted/stability path can have an
    // isolated invalid midpoint even though valid states exist on both sides.
    // A local scan preserves the physical unstable/stable bracket without
    // converting one solver miss into a missing boundary.
    double low = selectedLow;
    double high = selectedHigh;
    for (int level = 0; level < refinementLevels; ++level)
    {
        bool localHavePrevious = false;
        double localPreviousPressure = 0.0;
        OwStabilityState localPrevious;
        double refinedLow = std::numeric_limits<double>::quiet_NaN();
        double refinedHigh = std::numeric_limits<double>::quiet_NaN();

        for (int segment = 0; segment <= refinementSegments; ++segment)
        {
            const double fraction =
                static_cast<double>(segment) /
                static_cast<double>(refinementSegments);
            const double pressure = low + fraction * (high - low);
            const auto state =
                evaluateOwStability(flash, pressure, temperatureK, z);
            ++result.evaluations;
            if (!state.valid)
                continue;

            if (localHavePrevious &&
                localPrevious.gasUnstable && !state.gasUnstable)
            {
                // Keep the highest-pressure crossing within this bracket too.
                refinedLow = localPreviousPressure;
                refinedHigh = pressure;
            }
            localPreviousPressure = pressure;
            localPrevious = state;
            localHavePrevious = true;
        }

        if (!std::isfinite(refinedLow) || !std::isfinite(refinedHigh))
            break;
        low = refinedLow;
        high = refinedHigh;
    }

    result.found = std::isfinite(low) && std::isfinite(high) && high > low;
    if (result.found)
    {
        result.bracketLowMPa = low;
        result.bracketHighMPa = high;
        result.pressureMPa = 0.5 * (low + high);
    }
    return result;
}

struct BoundaryMetrics
{
    std::string convention;
    double waterMassFraction{};
    double waterMoleFraction{};
    std::size_t found{};
    double aadMPa{std::numeric_limits<double>::quiet_NaN()};
    double maxAbsErrorMPa{std::numeric_limits<double>::quiet_NaN()};
    double publishedAadDifferenceMPa{std::numeric_limits<double>::quiet_NaN()};
};

BoundaryMetrics evaluateBoundaryConvention(
    const Flash &flash,
    const std::vector<BoundaryPoint> &points,
    const std::string &convention,
    double waterMassFraction,
    std::ofstream &rows)
{
    constexpr double publishedAadMPa = 0.771;
    const Composition z = feedFromWaterMassFraction(waterMassFraction);

    BoundaryMetrics metrics;
    metrics.convention = convention;
    metrics.waterMassFraction = waterMassFraction;
    metrics.waterMoleFraction = z[water];

    double sumAbs = 0.0;
    double maxAbs = 0.0;
    for (const auto &point : points)
    {
        const auto predicted =
            findWlvWlBoundary(flash, point.temperatureK, z);
        const double error = predicted.found ?
            std::abs(predicted.pressureMPa - point.pressureMPa) :
            std::numeric_limits<double>::quiet_NaN();
        if (predicted.found)
        {
            ++metrics.found;
            sumAbs += error;
            maxAbs = std::max(maxAbs, error);
        }
        rows << convention << ',' << waterMassFraction << ',' << z[water]
             << ',' << point.temperatureK << ',' << point.pressureMPa
             << ',' << point.temperatureUncertaintyK << ','
             << point.pressureUncertaintyMPa << ','
             << (predicted.found ? 1 : 0) << ','
             << predicted.pressureMPa << ','
             << error << ',' << predicted.bracketLowMPa << ','
             << predicted.bracketHighMPa << ','
             << predicted.transitionCount << ','
             << predicted.evaluations << '\n';
    }

    if (metrics.found == points.size() && !points.empty())
    {
        metrics.aadMPa =
            sumAbs / static_cast<double>(points.size());
        metrics.maxAbsErrorMPa = maxAbs;
        metrics.publishedAadDifferenceMPa =
            std::abs(metrics.aadMPa - publishedAadMPa);
    }
    return metrics;
}

int run(
    const std::filesystem::path &referencePath,
    const std::filesystem::path &boundaryPath,
    const std::filesystem::path &outputDir)
{
    const auto points = readReference(referencePath);
    if (points.empty())
        throw std::runtime_error("Athabasca reference table contains no rows.");
    const auto boundaryPoints = readBoundaryReference(boundaryPath);
    if (boundaryPoints.empty())
        throw std::runtime_error(
            "Athabasca Figure-7 boundary table contains no rows.");
    std::filesystem::create_directories(outputDir);

    Eos eos = makeJiaCase1Cpa();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = static_cast<int>(water);
    options.maximumIterations = 240;
    options.maximumStabilityIterations = 160;
    Flash flash(eos, options);
    const Composition z = experimentalFeed();

    std::ofstream rows(outputDir / "water_solubility_comparison.csv");
    std::ofstream metrics(outputDir / "metrics.csv");
    std::ofstream gate(outputDir / "benchmark_gate.csv");
    std::ofstream boundaryRows(
        outputDir / "figure7_boundary_comparison.csv");
    std::ofstream boundaryMetrics(
        outputDir / "figure7_boundary_metrics.csv");
    if (!rows || !metrics || !gate || !boundaryRows || !boundaryMetrics)
        throw std::runtime_error("Cannot create Athabasca benchmark outputs.");

    boundaryRows << std::scientific << std::setprecision(12)
        << "feed_convention,water_mass_fraction,water_mole_fraction,"
           "T_K,experimental_P_MPa,T_uncertainty_K,P_uncertainty_MPa,"
           "boundary_found,predicted_P_MPa,abs_error_MPa,"
           "bracket_low_MPa,bracket_high_MPa,transition_count,evaluations\n";

    rows << std::scientific << std::setprecision(12);
    rows << "T_K,P_MPa,zH2O_feed,converged,phase_code,phase_count,"
            "beta_oil,beta_gas,beta_water,xH2O_oil,xH2O_water,"
            "experimental_xH2O_oil,published_Jia_CPA_xH2O_oil,"
            "ow_branch_abs_error_vs_experiment,"
            "ow_branch_abs_error_vs_published_CPA,"
            "max_material_closure,stability_valid,stability_stable,"
            "restricted_ow_converged,restricted_ow_phase_code,"
            "restricted_ow_phase_count,restricted_ow_xH2O_oil,"
            "restricted_ow_xH2O_water,restricted_ow_material_closure,"
            "restricted_ow_stability_valid,restricted_ow_stability_stable\n";

    std::size_t unrestrictedFailures = 0;
    std::size_t owBranchFailures = 0;
    std::size_t unrestrictedFailuresWithOwBranch = 0;
    double sumAbsExperimental = 0.0;
    double sumAbsPublished = 0.0;
    double maxAbsExperimental = 0.0;
    double maxUnrestrictedClosure = 0.0;
    double maxOwBranchClosure = 0.0;

    for (const auto &point : points)
    {
        const double pressure = point.pressureMPa * 1.0e6;

        // Unrestricted equilibrium is retained as a phase-topology diagnostic.
        // Table 5 is reported at WLV-WL transition points, so an exact
        // experimental T/P may fall on either side of the model-predicted
        // boundary when the predicted transition pressure is offset.
        const auto result = flash.flash(pressure, point.temperatureK, z);
        const bool hasOil = result.converged &&
            result.presence.contains(MPMC::CompositionalPhase::Oil);
        const bool hasWater = result.converged &&
            result.presence.contains(MPMC::CompositionalPhase::Water);
        double closure = std::numeric_limits<double>::quiet_NaN();
        bool stabilityValid = false;
        bool stabilityStable = false;
        double xOil = std::numeric_limits<double>::quiet_NaN();
        double xWater = std::numeric_limits<double>::quiet_NaN();

        if (result.converged)
        {
            closure = maxMaterialClosure(z, result);
            const auto stability = flash.stabilityTest(
                pressure, point.temperatureK, z,
                result.presence, result.composition);
            stabilityValid = stability.valid;
            stabilityStable = stability.stable;
            maxUnrestrictedClosure = std::max(
                maxUnrestrictedClosure, closure);
            if (hasOil)
                xOil = result.composition[
                    MPMC::phaseIndex(MPMC::CompositionalPhase::Oil)][water];
            if (hasWater)
                xWater = result.composition[
                    MPMC::phaseIndex(MPMC::CompositionalPhase::Water)][water];
        }

        const bool unrestrictedPass =
            result.converged && hasOil && hasWater &&
            stabilityValid && stabilityStable &&
            std::isfinite(closure) && closure <= 1.0e-8 &&
            std::isfinite(xOil);
        unrestrictedFailures += unrestrictedPass ? 0u : 1u;

        // Jia & Okuno Table 5 reports the water content of the bitumen-rich
        // liquid at experimental WLV-WL transition points.  Reproduce that
        // liquid branch consistently with an O+W restricted solve; global
        // stability is reported separately because an incipient vapor is the
        // defining missing phase at this boundary.
        const auto oilWater = MPMC::PhasePresence(
            static_cast<std::uint8_t>(
                MPMC::PhasePresence::oilBit |
                MPMC::PhasePresence::waterBit));
        const auto restricted = flash.flashRestricted(
            pressure, point.temperatureK, z, oilWater);

        const bool restrictedOwConverged = restricted.converged;
        const int restrictedOwPhaseCode = restricted.converged ?
            static_cast<int>(restricted.presence.bits()) : 0;
        const int restrictedOwPhaseCount = restricted.converged ?
            restricted.presence.count() : 0;
        double restrictedOwXOil = std::numeric_limits<double>::quiet_NaN();
        double restrictedOwXWater = std::numeric_limits<double>::quiet_NaN();
        double restrictedOwClosure = std::numeric_limits<double>::quiet_NaN();
        bool restrictedOwStabilityValid = false;
        bool restrictedOwStabilityStable = false;

        if (restricted.converged)
        {
            restrictedOwClosure = maxMaterialClosure(z, restricted);
            maxOwBranchClosure = std::max(
                maxOwBranchClosure, restrictedOwClosure);
            if (restricted.presence.contains(
                    MPMC::CompositionalPhase::Oil))
            {
                restrictedOwXOil = restricted.composition[
                    MPMC::phaseIndex(MPMC::CompositionalPhase::Oil)][water];
            }
            if (restricted.presence.contains(
                    MPMC::CompositionalPhase::Water))
            {
                restrictedOwXWater = restricted.composition[
                    MPMC::phaseIndex(MPMC::CompositionalPhase::Water)][water];
            }
            const auto restrictedStability = flash.stabilityTest(
                pressure, point.temperatureK, z,
                restricted.presence, restricted.composition);
            restrictedOwStabilityValid = restrictedStability.valid;
            restrictedOwStabilityStable = restrictedStability.stable;
        }

        const bool owBranchPass =
            restrictedOwConverged &&
            restrictedOwPhaseCode == static_cast<int>(oilWater.bits()) &&
            std::isfinite(restrictedOwClosure) &&
            restrictedOwClosure <= 1.0e-8 &&
            std::isfinite(restrictedOwXOil) &&
            std::isfinite(restrictedOwXWater);
        owBranchFailures += owBranchPass ? 0u : 1u;
        if (!unrestrictedPass && owBranchPass)
            ++unrestrictedFailuresWithOwBranch;

        const double errorExperimental = owBranchPass ?
            std::abs(restrictedOwXOil -
                     point.experimentalWaterInBitumen) :
            std::numeric_limits<double>::quiet_NaN();
        const double errorPublished = owBranchPass ?
            std::abs(restrictedOwXOil -
                     point.publishedCpaWaterInBitumen) :
            std::numeric_limits<double>::quiet_NaN();
        if (owBranchPass)
        {
            sumAbsExperimental += errorExperimental;
            sumAbsPublished += errorPublished;
            maxAbsExperimental = std::max(
                maxAbsExperimental, errorExperimental);
        }

        rows << point.temperatureK << ',' << point.pressureMPa << ','
             << z[water] << ',' << (result.converged ? 1 : 0) << ','
             << (result.converged ? static_cast<int>(result.presence.bits()) : 0)
             << ',' << (result.converged ? result.presence.count() : 0) << ','
             << result.phaseMoleFraction[0] << ','
             << result.phaseMoleFraction[1] << ','
             << result.phaseMoleFraction[2] << ','
             << xOil << ',' << xWater << ','
             << point.experimentalWaterInBitumen << ','
             << point.publishedCpaWaterInBitumen << ','
             << errorExperimental << ',' << errorPublished << ','
             << closure << ',' << (stabilityValid ? 1 : 0) << ','
             << (stabilityStable ? 1 : 0) << ','
             << (restrictedOwConverged ? 1 : 0) << ','
             << restrictedOwPhaseCode << ',' << restrictedOwPhaseCount << ','
             << restrictedOwXOil << ',' << restrictedOwXWater << ','
             << restrictedOwClosure << ','
             << (restrictedOwStabilityValid ? 1 : 0) << ','
             << (restrictedOwStabilityStable ? 1 : 0) << '\n';
    }

    const std::size_t unrestrictedPassed =
        points.size() - unrestrictedFailures;
    const std::size_t owBranchPassed =
        points.size() - owBranchFailures;
    const double maeExperimental = owBranchPassed > 0 ?
        sumAbsExperimental / static_cast<double>(owBranchPassed) :
        std::numeric_limits<double>::quiet_NaN();
    const double maePublished = owBranchPassed > 0 ?
        sumAbsPublished / static_cast<double>(owBranchPassed) :
        std::numeric_limits<double>::quiet_NaN();

    metrics << std::scientific << std::setprecision(12);
    metrics << "reference_points,unrestricted_pass,unrestricted_fail,"
               "ow_branch_pass,ow_branch_fail,"
               "unrestricted_failures_with_ow_branch,"
               "mae_vs_experiment,mae_vs_published_Jia_CPA,"
               "max_abs_error_vs_experiment,max_unrestricted_material_closure,"
               "max_ow_branch_material_closure\n";
    metrics << points.size() << ',' << unrestrictedPassed << ','
            << unrestrictedFailures << ',' << owBranchPassed << ','
            << owBranchFailures << ','
            << unrestrictedFailuresWithOwBranch << ','
            << maeExperimental << ',' << maePublished << ','
            << maxAbsExperimental << ',' << maxUnrestrictedClosure << ','
            << maxOwBranchClosure << '\n';

    const auto jiaCaptionBoundary = evaluateBoundaryConvention(
        flash, boundaryPoints, "JIA_FIGURE7_CAPTION",
        0.559, boundaryRows);
    const auto amaniTableBoundary = evaluateBoundaryConvention(
        flash, boundaryPoints, "AMANI_TABLE_2_5",
        0.441, boundaryRows);

    constexpr double publishedBoundaryAadMPa = 0.771;
    // Jia & Okuno report only the aggregate boundary AAD, not the individual
    // calculated Figure-7 pressures.  Until the exact source composition and
    // CPA convention are reconciled, the published AAD is an audit target,
    // not a hard acceptance threshold.  Do not tune parameters to reproduce
    // the magnitude of an aggregate error statistic.

    const auto writeBoundaryMetric = [&](const BoundaryMetrics &m) {
        boundaryMetrics << std::scientific << std::setprecision(12)
            << m.convention << ',' << m.waterMassFraction << ','
            << m.waterMoleFraction << ',' << m.found << ','
            << boundaryPoints.size() << ',' << m.aadMPa << ','
            << m.maxAbsErrorMPa << ',' << publishedBoundaryAadMPa << ','
            << m.publishedAadDifferenceMPa << '\n';
    };
    boundaryMetrics
        << "feed_convention,water_mass_fraction,water_mole_fraction,"
           "boundary_points_found,boundary_points_total,aad_MPa,"
           "max_abs_error_MPa,published_Jia_aad_MPa,"
           "abs_difference_from_published_aad_MPa\n";
    writeBoundaryMetric(jiaCaptionBoundary);
    writeBoundaryMetric(amaniTableBoundary);

    const bool boundaryRowsComplete =
        jiaCaptionBoundary.found == boundaryPoints.size() &&
        amaniTableBoundary.found == boundaryPoints.size();

    // Table 5 is a branch-composition benchmark at experimental WLV-WL
    // transition points.  Therefore Stage 1 hard-gates reproducibility of the
    // O+W liquid branch, while unrestricted phase topology remains diagnostic
    // until the separate Figure-7 boundary benchmark is implemented.
    const bool owBranchGate = owBranchFailures == 0;
    gate << "gate,status,value,criterion\n";
    gate << "CPA_ATHABASCA_TABLE5_OW_BRANCH,"
         << (owBranchGate ? "PASS" : "FAIL") << ','
         << owBranchPassed << '/' << points.size()
         << ",all Table-5 T/P rows must reproduce a finite O+W branch with closure<=1e-8\n";
    gate << "CPA_ATHABASCA_UNRESTRICTED_EQUILIBRIUM,OBSERVE,"
         << unrestrictedPassed << '/' << points.size()
         << ",diagnostic only because Table 5 rows are WLV-WL transition points; validate topology against Figure 7 separately\n";
    gate << "CPA_ATHABASCA_PARITY,OBSERVE,"
         << std::scientific << std::setprecision(12) << maeExperimental
         << ",O+W branch MAE versus experiment; report-only until a preregistered parity tolerance is fixed\n";
    gate << "CPA_ATHABASCA_FIGURE7_BOUNDARY_ROWS,"
         << (boundaryRowsComplete ? "PASS" : "FAIL") << ','
         << (jiaCaptionBoundary.found + amaniTableBoundary.found) << '/'
         << (2 * boundaryPoints.size())
         << ",both conflicting source feed conventions must yield one traceable WLV-WL boundary at every experimental temperature\n";
    gate << "CPA_ATHABASCA_FIGURE7_JIA_CAPTION_AAD,OBSERVE,"
         << jiaCaptionBoundary.aadMPa
         << ",audit against published 0.771 MPa AAD; Jia caption composition conflicts with the original Amani table\n";
    gate << "CPA_ATHABASCA_FIGURE7_AMANI_TABLE_AAD,OBSERVE,"
         << amaniTableBoundary.aadMPa
         << ",audit against published 0.771 MPa AAD using the original Amani Table-2.5 composition label\n";
    gate << "CPA_ATHABASCA_FIGURE7_PUBLISHED_AAD_DIFFERENCE,OBSERVE,"
         << std::min(
                jiaCaptionBoundary.publishedAadDifferenceMPa,
                amaniTableBoundary.publishedAadDifferenceMPa)
         << ",no hard tolerance until the source composition and CPA convention are reconciled; never tune parameters to match aggregate AAD\n";

    std::cout << "Athabasca CPA proxy: O+W branch "
              << owBranchPassed << '/' << points.size()
              << ", unrestricted equilibrium " << unrestrictedPassed << '/'
              << points.size() << ", Table5 MAE(exp)=" << std::scientific
              << maeExperimental << ", Figure7 AAD Jia-caption="
              << jiaCaptionBoundary.aadMPa << " MPa, Amani-table="
              << amaniTableBoundary.aadMPa
              << " MPa, published=0.771 MPa\n";
    return (owBranchGate && boundaryRowsComplete) ? 0 : 2;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 4)
        {
            std::cerr
                << "usage: cpa_athabasca_bitumen_water "
                   "SOLUBILITY_CSV BOUNDARY_CSV OUTPUT_DIR\n";
            return 1;
        }
        return run(argv[1], argv[2], argv[3]);
    }
    catch (const std::exception &error)
    {
        std::cerr << "Athabasca CPA proxy benchmark failed: "
                  << error.what() << '\n';
        return 1;
    }
}
