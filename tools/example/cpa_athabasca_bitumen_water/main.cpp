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

Composition experimentalFeed()
{
    // Figure 7 of Jia & Okuno Case 1: 55.9 wt% water + 44.1 wt%
    // Athabasca bitumen.  Table 4 z-values are normalized here because the
    // published rounded values sum to 0.999.
    const Composition oil = bitumenComposition();
    constexpr std::array<double, 5> mwKgMol{
        0.01801528, 0.35243, 0.53903, 0.70684, 0.91619};
    double oilMw = 0.0;
    for (std::size_t i = 1; i < oil.size(); ++i)
        oilMw += oil[i] * mwKgMol[i];

    const double nWater = 0.559 / mwKgMol[0];
    const double nOil = 0.441 / oilMw;
    const double xWater = nWater / (nWater + nOil);

    Composition z{};
    z[0] = xWater;
    for (std::size_t i = 1; i < z.size(); ++i)
        z[i] = (1.0 - xWater) * oil[i];
    return z;
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

    // Standard 4C water.
    cpa.a0[water] = MPMC::StandardCpaWater4C::a0;
    cpa.b[water] = MPMC::StandardCpaWater4C::b;
    cpa.c1[water] = MPMC::StandardCpaWater4C::c1;
    cpa.associationEnergy[water] = MPMC::StandardCpaWater4C::epsilon;
    cpa.associationVolume[water] = MPMC::StandardCpaWater4C::beta;
    cpa.donorSites[water] = MPMC::StandardCpaWater4C::donorSites;
    cpa.acceptorSites[water] = MPMC::StandardCpaWater4C::acceptorSites;

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
            0.5 * MPMC::StandardCpaWater4C::epsilon;
        cpa.crossAssociationVolume[water][pc] =
            std::sqrt(MPMC::StandardCpaWater4C::beta * 0.07);
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

int run(
    const std::filesystem::path &referencePath,
    const std::filesystem::path &outputDir)
{
    const auto points = readReference(referencePath);
    if (points.empty())
        throw std::runtime_error("Athabasca reference table contains no rows.");
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
    if (!rows || !metrics || !gate)
        throw std::runtime_error("Cannot create Athabasca benchmark outputs.");

    rows << std::scientific << std::setprecision(12);
    rows << "T_K,P_MPa,zH2O_feed,converged,phase_code,phase_count,"
            "beta_oil,beta_gas,beta_water,xH2O_oil,xH2O_water,"
            "experimental_xH2O_oil,published_Jia_CPA_xH2O_oil,"
            "abs_error_vs_experiment,abs_error_vs_published_CPA,"
            "max_material_closure,stability_valid,stability_stable,"
            "restricted_ow_converged,restricted_ow_phase_code,"
            "restricted_ow_phase_count,restricted_ow_xH2O_oil,"
            "restricted_ow_xH2O_water,restricted_ow_material_closure,"
            "restricted_ow_stability_valid,restricted_ow_stability_stable\n";

    std::size_t structuralFailures = 0;
    std::size_t restrictedOwRecoveries = 0;
    double sumAbsExperimental = 0.0;
    double sumAbsPublished = 0.0;
    double maxAbsExperimental = 0.0;
    double maxClosure = 0.0;

    for (const auto &point : points)
    {
        const double pressure = point.pressureMPa * 1.0e6;
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

        bool restrictedOwConverged = false;
        int restrictedOwPhaseCode = 0;
        int restrictedOwPhaseCount = 0;
        double restrictedOwXOil = std::numeric_limits<double>::quiet_NaN();
        double restrictedOwXWater = std::numeric_limits<double>::quiet_NaN();
        double restrictedOwClosure = std::numeric_limits<double>::quiet_NaN();
        bool restrictedOwStabilityValid = false;
        bool restrictedOwStabilityStable = false;

        if (result.converged)
        {
            closure = maxMaterialClosure(z, result);
            const auto stability = flash.stabilityTest(
                pressure, point.temperatureK, z,
                result.presence, result.composition);
            stabilityValid = stability.valid;
            stabilityStable = stability.stable;
            maxClosure = std::max(maxClosure, closure);
            if (hasOil)
                xOil = result.composition[
                    MPMC::phaseIndex(MPMC::CompositionalPhase::Oil)][water];
            if (hasWater)
                xWater = result.composition[
                    MPMC::phaseIndex(MPMC::CompositionalPhase::Water)][water];
        }
        else
        {
            // Diagnostic only: if the unrestricted active-set path fails,
            // ask whether a physically certified Oil+Water solution exists.
            // This does not change the structural gate and does not tune any
            // CPA parameter.  A stable restricted solution identifies a
            // solver-path robustness issue rather than missing equilibrium.
            const auto oilWater = MPMC::PhasePresence(
                static_cast<std::uint8_t>(
                    MPMC::PhasePresence::oilBit |
                    MPMC::PhasePresence::waterBit));
            const auto restricted = flash.flashRestricted(
                pressure, point.temperatureK, z, oilWater);
            restrictedOwConverged = restricted.converged;
            if (restricted.converged)
            {
                restrictedOwPhaseCode =
                    static_cast<int>(restricted.presence.bits());
                restrictedOwPhaseCount = restricted.presence.count();
                restrictedOwClosure = maxMaterialClosure(z, restricted);
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
                if (restrictedOwPhaseCode ==
                        static_cast<int>(oilWater.bits()) &&
                    restrictedOwStabilityValid &&
                    restrictedOwStabilityStable &&
                    std::isfinite(restrictedOwClosure) &&
                    restrictedOwClosure <= 1.0e-8)
                {
                    ++restrictedOwRecoveries;
                }
            }
        }

        const bool structuralPass =
            result.converged && hasOil && hasWater &&
            stabilityValid && stabilityStable &&
            std::isfinite(closure) && closure <= 1.0e-8 &&
            std::isfinite(xOil);
        structuralFailures += structuralPass ? 0u : 1u;

        const double errorExperimental = structuralPass ?
            std::abs(xOil - point.experimentalWaterInBitumen) :
            std::numeric_limits<double>::quiet_NaN();
        const double errorPublished = structuralPass ?
            std::abs(xOil - point.publishedCpaWaterInBitumen) :
            std::numeric_limits<double>::quiet_NaN();
        if (structuralPass)
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

    const std::size_t passed = points.size() - structuralFailures;
    const double maeExperimental = passed > 0 ?
        sumAbsExperimental / static_cast<double>(passed) :
        std::numeric_limits<double>::quiet_NaN();
    const double maePublished = passed > 0 ?
        sumAbsPublished / static_cast<double>(passed) :
        std::numeric_limits<double>::quiet_NaN();

    metrics << std::scientific << std::setprecision(12);
    metrics << "reference_points,structural_pass,structural_fail,"
               "restricted_ow_recoveries,mae_vs_experiment,"
               "mae_vs_published_Jia_CPA,max_abs_error_vs_experiment,"
               "max_material_closure\n";
    metrics << points.size() << ',' << passed << ',' << structuralFailures
            << ',' << restrictedOwRecoveries << ','
            << maeExperimental << ',' << maePublished << ','
            << maxAbsExperimental << ',' << maxClosure << '\n';

    // Stage 1 is intentionally structural.  A parity threshold is not
    // invented before the production implementation has been compared once
    // against the published Case-1 convention.
    const bool structuralGate = structuralFailures == 0;
    gate << "gate,status,value,criterion\n";
    gate << "CPA_ATHABASCA_STRUCTURAL,"
         << (structuralGate ? "PASS" : "FAIL") << ','
         << passed << '/' << points.size()
         << ",all rows converge with stable Oil+Water and closure<=1e-8\n";
    gate << "CPA_ATHABASCA_PARITY,OBSERVE,"
         << std::scientific << std::setprecision(12) << maeExperimental
         << ",report-only until convention audit fixes a preregistered tolerance\n";
    gate << "CPA_ATHABASCA_RESTRICTED_OW_DIAGNOSTIC,OBSERVE,"
         << restrictedOwRecoveries
         << ",count unrestricted failures that have a stable certified O+W restricted solution\n";

    std::cout << "Athabasca CPA proxy: structural " << passed << '/'
              << points.size() << ", restricted O+W recoveries="
              << restrictedOwRecoveries << ", MAE(exp)=" << std::scientific
              << maeExperimental << ", MAE(Jia CPA)=" << maePublished
              << ", max closure=" << maxClosure << '\n';
    return structuralGate ? 0 : 2;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            std::cerr
                << "usage: cpa_athabasca_bitumen_water REFERENCE_CSV OUTPUT_DIR\n";
            return 1;
        }
        return run(argv[1], argv[2]);
    }
    catch (const std::exception &error)
    {
        std::cerr << "Athabasca CPA proxy benchmark failed: "
                  << error.what() << '\n';
        return 1;
    }
}
