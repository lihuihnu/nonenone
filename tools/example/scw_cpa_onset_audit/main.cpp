/**
 * @file main.cpp
 * @brief Narrow CPA phase-onset continuation audit with frozen parameters.
 *
 * This diagnostic never changes CPA a0/b/c1, association or BIP parameters.
 * It compares the production unrestricted flash against restricted O+W
 * equilibria, stability certificates and fine seeded continuation paths.
 */
#include <common/math.hpp>
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
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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

constexpr std::size_t N = 5;
constexpr std::size_t water = 0;
constexpr double temperature = 653.15;
constexpr double failureWaterFraction = 0.7610625;
constexpr double missing = std::numeric_limits<double>::quiet_NaN();

inline const std::array<std::string, N> names{
    "H2O", "OIL_GASOLINE", "OIL_DIESEL", "OIL_MIDDLE", "OIL_HEAVY"};

struct Csv
{
    std::map<std::string, std::size_t> column;
    std::vector<std::vector<std::string>> row;
};

struct Parameters
{
    std::array<double, N> tc{};
    std::array<double, N> pc{};
    std::array<double, N> vc{};
    std::array<double, N> omega{};
    std::array<double, N> mw{};
    std::array<double, N> a0{};
    std::array<double, N> b{};
    std::array<double, N> c1{};
    std::array<double, N> associationEnergy{};
    std::array<double, N> associationVolume{};
    std::array<int, N> donorSites{};
    std::array<int, N> acceptorSites{};
    Matrix kij{};
    std::array<double, 4> baseOilRatio{};
};

struct Candidate
{
    bool converged{false};
    int phaseCode{0};
    int phaseCount{0};
    int iterations{0};
    bool stabilityValid{false};
    bool stabilityStable{false};
    bool oilMissingUnstable{false};
    bool gasMissingUnstable{false};
    bool waterMissingUnstable{false};
    double trialOil{missing};
    double trialGas{missing};
    double trialWater{missing};
    double betaOil{missing};
    double betaWater{missing};
    double xOilWater{missing};
    double xWaterWater{missing};
    bool roleConsistent{false};
    double massClosure{missing};
    double logFugacitySpread{missing};
    Flash::Result result{};
};

struct RootProbe
{
    bool oilOk{false};
    bool gasOk{false};
    bool waterOk{false};
    double oilZ{missing};
    double gasZ{missing};
    double waterZ{missing};
    MPMC::CubicEquationOfState<Indices>::ThermodynamicProfile profile{};
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

Csv readCsv(const std::filesystem::path &path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open CSV: " + path.string());
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("Empty CSV: " + path.string());
    Csv result;
    const auto header = splitCsv(line);
    for (std::size_t i = 0; i < header.size(); ++i)
        result.column.emplace(header[i], i);
    while (std::getline(in, line))
        if (!line.empty())
            result.row.push_back(splitCsv(line));
    return result;
}

const std::string &field(
    const Csv &csv,
    const std::vector<std::string> &row,
    const std::string &name)
{
    return row.at(csv.column.at(name));
}

std::size_t componentIndex(const std::string &name)
{
    for (std::size_t i = 0; i < names.size(); ++i)
        if (names[i] == name)
            return i;
    throw std::runtime_error("Unknown component: " + name);
}

void normalize(Composition &x)
{
    double sum = 0.0;
    for (double value : x)
    {
        if (!std::isfinite(value) || value < 0.0)
            throw std::runtime_error("Invalid composition.");
        sum += value;
    }
    if (!(sum > 0.0))
        throw std::runtime_error("Empty composition.");
    for (double &value : x)
        value /= sum;
}

Parameters readParameters(const std::filesystem::path &caseDir)
{
    Parameters p;
    const auto pure = readCsv(
        caseDir / "cpa_parameters/cpa_pure_parameters_380c.csv");
    for (const auto &row : pure.row)
    {
        const std::size_t i = componentIndex(field(pure, row, "component"));
        p.tc[i] = std::stod(field(pure, row, "Tc_K"));
        p.pc[i] = std::stod(field(pure, row, "Pc_MPa")) * 1.0e6;
        p.vc[i] = std::stod(field(pure, row, "Vc_cm3_mol")) * 1.0e-6;
        p.omega[i] = std::stod(field(pure, row, "omega"));
        p.mw[i] = std::stod(field(pure, row, "MW_g_mol")) * 1.0e-3;
        p.a0[i] = std::stod(field(pure, row, "cpa_a0_Pa_m6_mol2"));
        p.b[i] = std::stod(field(pure, row, "cpa_b_m3_mol"));
        p.c1[i] = std::stod(field(pure, row, "cpa_c1"));
    }

    const auto association = readCsv(
        caseDir / "cpa_parameters/cpa_association_scheme.csv");
    for (const auto &row : association.row)
    {
        const std::size_t i =
            componentIndex(field(association, row, "component"));
        p.donorSites[i] =
            std::stoi(field(association, row, "donor_sites"));
        p.acceptorSites[i] =
            std::stoi(field(association, row, "acceptor_sites"));
        p.associationEnergy[i] =
            std::stod(field(association, row, "association_energy_J_mol"));
        p.associationVolume[i] =
            std::stod(field(association, row, "association_volume"));
    }

    const auto bip = readCsv(
        caseDir / "cpa_parameters/cpa_binary_matrix_screening.csv");
    for (const auto &row : bip.row)
    {
        const std::size_t i = componentIndex(field(bip, row, "component_i"));
        const std::size_t j = componentIndex(field(bip, row, "component_j"));
        const double value =
            std::stod(field(bip, row, "kij_constant_initial"));
        p.kij[i][j] = p.kij[j][i] = value;
    }

    const auto scan = readCsv(caseDir / "pvt_acceptance/composition_scan.csv");
    bool found = false;
    for (const auto &row : scan.row)
    {
        if (field(scan, row, "family") != "BASE")
            continue;
        const double zw = std::stod(field(scan, row, "z_H2O"));
        if (std::abs(zw - 0.20) > 1.0e-10)
            continue;
        const double oil = 1.0 - zw;
        p.baseOilRatio = {
            std::stod(field(scan, row, "z_OIL_GASOLINE")) / oil,
            std::stod(field(scan, row, "z_OIL_DIESEL")) / oil,
            std::stod(field(scan, row, "z_OIL_MIDDLE")) / oil,
            std::stod(field(scan, row, "z_OIL_HEAVY")) / oil};
        double sum = 0.0;
        for (double value : p.baseOilRatio)
            sum += value;
        for (double &value : p.baseOilRatio)
            value /= sum;
        found = true;
        break;
    }
    if (!found)
        throw std::runtime_error("Cannot recover BASE oil ratio.");
    return p;
}

Composition compositionAt(
    const Parameters &p,
    double waterFraction)
{
    Composition z{};
    z[0] = waterFraction;
    const double oil = 1.0 - waterFraction;
    for (std::size_t i = 0; i < 4; ++i)
        z[i + 1] = oil * p.baseOilRatio[i];
    normalize(z);
    return z;
}

Eos makeCpa(const Parameters &p)
{
    Eos eos(
        0.42748, 0.08664,
        Mixture(p.tc, p.pc, p.vc, p.omega, p.mw, p.kij),
        1, 1.0, 0.0, 1.0e-30);
    Eos::CubicPlusAssociationOptions cpa;
    cpa.a0 = p.a0;
    cpa.b = p.b;
    cpa.c1 = p.c1;
    cpa.associationEnergy = p.associationEnergy;
    cpa.associationVolume = p.associationVolume;
    cpa.donorSites = p.donorSites;
    cpa.acceptorSites = p.acceptorSites;
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::SoaveRedlichKwong;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    eos.configureCpaTemperatureCache(temperature);
    eos.setThermodynamicProfilerEnabled(true);
    return eos;
}

double closure(const Composition &z, const Flash::Result &result)
{
    double maximum = 0.0;
    for (std::size_t c = 0; c < N; ++c)
    {
        double reconstructed = 0.0;
        for (std::size_t phase = 0; phase < 3; ++phase)
            reconstructed += result.phaseMoleFraction[phase] *
                result.composition[phase][c];
        maximum = std::max(maximum, std::abs(reconstructed - z[c]));
    }
    return maximum;
}

double fugacitySpread(
    const Eos &eos,
    double pressure,
    const Flash::Result &result)
{
    if (result.presence.count() <= 1)
        return 0.0;
    double maximum = 0.0;
    for (std::size_t c = 0; c < N; ++c)
    {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        int count = 0;
        for (std::size_t phase = 0; phase < 3; ++phase)
        {
            const auto role =
                static_cast<MPMC::CompositionalPhase>(phase);
            if (!result.presence.contains(role))
                continue;
            const auto thermo = eos.phaseResult(
                pressure, temperature, result.composition[phase], role, false);
            const double f = std::max(thermo.fugacity[c], 1.0e-300);
            lo = std::min(lo, std::log(f));
            hi = std::max(hi, std::log(f));
            ++count;
        }
        if (count > 1)
            maximum = std::max(maximum, hi - lo);
    }
    return maximum;
}

Candidate certify(
    const Eos &eos,
    const Flash &flash,
    double pressure,
    const Composition &z,
    Flash::Result result)
{
    Candidate candidate;
    candidate.result = result;
    candidate.converged = result.converged;
    candidate.iterations = result.iterations;
    if (!result.converged)
        return candidate;

    candidate.phaseCode = static_cast<int>(result.presence.bits());
    candidate.phaseCount = result.presence.count();
    const auto stability = flash.stabilityTest(
        pressure, temperature, z, result.presence, result.composition);
    candidate.stabilityValid = stability.valid;
    candidate.stabilityStable = stability.stable;
    candidate.oilMissingUnstable =
        stability.missingPhaseUnstable[0];
    candidate.gasMissingUnstable =
        stability.missingPhaseUnstable[1];
    candidate.waterMissingUnstable =
        stability.missingPhaseUnstable[2];
    candidate.trialOil = stability.trialSum[0];
    candidate.trialGas = stability.trialSum[1];
    candidate.trialWater = stability.trialSum[2];
    candidate.betaOil = result.phaseMoleFraction[0];
    candidate.betaWater = result.phaseMoleFraction[2];
    candidate.xOilWater = result.composition[0][water];
    candidate.xWaterWater = result.composition[2][water];
    candidate.roleConsistent =
        !result.presence.contains(MPMC::CompositionalPhase::Water) ||
        !result.presence.contains(MPMC::CompositionalPhase::Oil) ||
        candidate.xWaterWater + 1.0e-10 >= candidate.xOilWater;
    candidate.massClosure = closure(z, result);
    candidate.logFugacitySpread = fugacitySpread(eos, pressure, result);
    return candidate;
}

Flash::Result waterOnlyResult(
    const Flash &flash,
    double pressure,
    const Composition &z)
{
    return flash.flashRestricted(
        pressure, temperature, z, MPMC::PhasePresence::waterOnly());
}

MPMC::PhasePresence oilWaterPresence()
{
    return MPMC::PhasePresence(static_cast<std::uint8_t>(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit));
}

RootProbe probeRoots(
    Eos &eos,
    double pressure,
    const Composition &z)
{
    RootProbe probe;
    eos.resetThermodynamicProfile();
    try
    {
        const auto phase = eos.phaseResult(
            pressure, temperature, z, MPMC::CompositionalPhase::Oil, false);
        probe.oilOk = std::isfinite(phase.compressibility) &&
            phase.compressibility > 0.0;
        probe.oilZ = phase.compressibility;
    }
    catch (const std::exception &) {}
    try
    {
        const auto phase = eos.phaseResult(
            pressure, temperature, z, MPMC::CompositionalPhase::Gas, false);
        probe.gasOk = std::isfinite(phase.compressibility) &&
            phase.compressibility > 0.0;
        probe.gasZ = phase.compressibility;
    }
    catch (const std::exception &) {}
    try
    {
        const auto phase = eos.phaseResult(
            pressure, temperature, z, MPMC::CompositionalPhase::Water, false);
        probe.waterOk = std::isfinite(phase.compressibility) &&
            phase.compressibility > 0.0;
        probe.waterZ = phase.compressibility;
    }
    catch (const std::exception &) {}
    probe.profile = eos.thermodynamicProfile();
    return probe;
}

void writeCandidate(
    std::ostream &out,
    const std::string &method,
    double waterFraction,
    double pressureMPa,
    const Candidate &c)
{
    out << method << ',' << waterFraction << ',' << pressureMPa << ','
        << (c.converged ? 1 : 0) << ',' << c.phaseCode << ','
        << c.phaseCount << ',' << c.iterations << ','
        << (c.stabilityValid ? 1 : 0) << ','
        << (c.stabilityStable ? 1 : 0) << ','
        << (c.oilMissingUnstable ? 1 : 0) << ','
        << (c.gasMissingUnstable ? 1 : 0) << ','
        << (c.waterMissingUnstable ? 1 : 0) << ','
        << c.trialOil << ',' << c.trialGas << ',' << c.trialWater << ','
        << c.betaOil << ',' << c.betaWater << ','
        << c.xOilWater << ',' << c.xWaterWater << ','
        << (c.roleConsistent ? 1 : 0) << ','
        << c.massClosure << ',' << c.logFugacitySpread << '\n';
}

double compositionDistance(
    const Flash::Result &a,
    const Flash::Result &b)
{
    double distance = 0.0;
    for (std::size_t phase : {std::size_t(0), std::size_t(2)})
        for (std::size_t c = 0; c < N; ++c)
            distance += std::abs(
                a.composition[phase][c] - b.composition[phase][c]);
    return distance;
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            std::cerr
                << "Usage: scw_cpa_onset_audit CASE_DIR OUTPUT_DIR\n";
            return 2;
        }
        const std::filesystem::path caseDir(argv[1]);
        const std::filesystem::path output(argv[2]);
        if (std::filesystem::exists(output))
            throw std::runtime_error("Output directory already exists.");
        std::filesystem::create_directories(output);

        const Parameters parameters = readParameters(caseDir);
        Eos eos = makeCpa(parameters);
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = static_cast<int>(water);
        options.maximumIterations = 180;
        options.maximumStabilityIterations = 120;
        const Flash flash(eos, options);

        std::ofstream diagnosis(output / "failure_state_diagnostics.csv");
        diagnosis << std::setprecision(17)
            << "method,z_H2O,pressure_MPa,converged,phase_code,phase_count,"
               "iterations,stability_valid,stability_stable,"
               "missing_oil_unstable,missing_gas_unstable,missing_water_unstable,"
               "trial_sum_oil,trial_sum_gas,trial_sum_water,beta_oil,beta_water,"
               "xH2O_oil,xH2O_water,role_consistent,mass_closure,"
               "max_log_fugacity_spread\n";

        std::ofstream roots(output / "cpa_root_association_probe.csv");
        roots << std::setprecision(17)
            << "z_H2O,pressure_MPa,oil_root_ok,gas_root_ok,water_root_ok,"
               "oil_Z,gas_Z,water_Z,cpa_phase_result_calls,cpa_density_root_calls,"
               "cpa_density_residual_evaluations,cpa_density_bisection_iterations,"
               "cpa_association_calls,cpa_association_analytic_calls,"
               "cpa_association_iterative_calls,cpa_association_iterations,"
               "water_only_analytic_fast_path\n";

        const std::array<double, 4> keyPressures{
            25.75, 26.00, 26.25, 26.50};
        const Composition failureZ =
            compositionAt(parameters, failureWaterFraction);

        for (double pressureMPa : keyPressures)
        {
            const double pressure = pressureMPa * 1.0e6;
            eos.resetThermodynamicProfile();
            writeCandidate(
                diagnosis, "unrestricted", failureWaterFraction, pressureMPa,
                certify(eos, flash, pressure, failureZ,
                        flash.flash(pressure, temperature, failureZ)));

            writeCandidate(
                diagnosis, "water_only", failureWaterFraction, pressureMPa,
                certify(eos, flash, pressure, failureZ,
                        waterOnlyResult(flash, pressure, failureZ)));

            writeCandidate(
                diagnosis, "oil_water_unseeded", failureWaterFraction, pressureMPa,
                certify(eos, flash, pressure, failureZ,
                        flash.flashRestricted(
                            pressure, temperature, failureZ, oilWaterPresence())));

            const auto probe = probeRoots(eos, pressure, failureZ);
            roots << failureWaterFraction << ',' << pressureMPa << ','
                << (probe.oilOk ? 1 : 0) << ','
                << (probe.gasOk ? 1 : 0) << ','
                << (probe.waterOk ? 1 : 0) << ','
                << probe.oilZ << ',' << probe.gasZ << ',' << probe.waterZ << ','
                << probe.profile.cpaPhaseResultCalls << ','
                << probe.profile.cpaDensityRootCalls << ','
                << probe.profile.cpaDensityResidualEvaluations << ','
                << probe.profile.cpaDensityBisectionIterations << ','
                << probe.profile.cpaAssociationCalls << ','
                << probe.profile.cpaAssociationAnalyticCalls << ','
                << probe.profile.cpaAssociationIterativeCalls << ','
                << probe.profile.cpaAssociationIterations << ','
                << (probe.profile.cpaWaterOnlyAnalyticFastPathAvailable ? 1 : 0)
                << '\n';
        }

        // Fine pressure continuation from the known O+W side.
        std::ofstream pressurePath(
            output / "seeded_pressure_continuation.csv");
        pressurePath << diagnosis.rdbuf();
        // Re-open with an explicit header because streambuf reuse does not
        // copy already-written bytes.
        pressurePath.close();
        pressurePath.open(
            output / "seeded_pressure_continuation.csv",
            std::ios::out | std::ios::trunc);
        pressurePath << std::setprecision(17)
            << "method,z_H2O,pressure_MPa,converged,phase_code,phase_count,"
               "iterations,stability_valid,stability_stable,"
               "missing_oil_unstable,missing_gas_unstable,missing_water_unstable,"
               "trial_sum_oil,trial_sum_gas,trial_sum_water,beta_oil,beta_water,"
               "xH2O_oil,xH2O_water,role_consistent,mass_closure,"
               "max_log_fugacity_spread\n";

        Flash::Result pressureSeed = flash.flashRestricted(
            27.0e6, temperature, failureZ, oilWaterPresence());
        if (!pressureSeed.converged)
            throw std::runtime_error(
                "CPA audit could not establish the 27 MPa O+W seed.");
        for (int step = 0; step <= 60; ++step)
        {
            const double pressureMPa = 27.0 - 0.025 * step;
            const double pressure = pressureMPa * 1.0e6;
            Flash::Result result = flash.flashRestricted(
                pressure, temperature, failureZ,
                oilWaterPresence(), pressureSeed.composition);
            const Candidate certified =
                certify(eos, flash, pressure, failureZ, result);
            writeCandidate(
                pressurePath, "seeded_high_to_low",
                failureWaterFraction, pressureMPa, certified);
            if (result.converged)
                pressureSeed = result;
        }

        // Fine composition continuation at the two failed pressures.
        std::ofstream compositionPath(
            output / "seeded_composition_continuation.csv");
        compositionPath << std::setprecision(17)
            << "method,z_H2O,pressure_MPa,converged,phase_code,phase_count,"
               "iterations,stability_valid,stability_stable,"
               "missing_oil_unstable,missing_gas_unstable,missing_water_unstable,"
               "trial_sum_oil,trial_sum_gas,trial_sum_water,beta_oil,beta_water,"
               "xH2O_oil,xH2O_water,role_consistent,mass_closure,"
               "max_log_fugacity_spread\n";

        for (double pressureMPa : {26.0, 26.25})
        {
            const double pressure = pressureMPa * 1.0e6;
            double zw = 0.775;
            Composition z = compositionAt(parameters, zw);
            Flash::Result seed = flash.flashRestricted(
                pressure, temperature, z, oilWaterPresence());
            if (!seed.converged)
                throw std::runtime_error(
                    "CPA audit could not establish the high-water O+W seed.");
            for (int step = 0; step <= 140; ++step)
            {
                zw = 0.775 - 0.00025 * step;
                z = compositionAt(parameters, zw);
                const Flash::Result result = flash.flashRestricted(
                    pressure, temperature, z,
                    oilWaterPresence(), seed.composition);
                const Candidate certified =
                    certify(eos, flash, pressure, z, result);
                writeCandidate(
                    compositionPath, "seeded_high_water_to_low",
                    zw, pressureMPa, certified);
                if (result.converged)
                    seed = result;
            }
        }

        // Role-canonicalization test using a converged O+W seed and a
        // deliberately swapped liquid-slot seed.
        std::ofstream role(output / "role_canonicalization_probe.csv");
        role << std::setprecision(17)
            << "pressure_MPa,normal_converged,swapped_converged,"
               "normal_xH2O_oil,normal_xH2O_water,"
               "swapped_xH2O_oil,swapped_xH2O_water,"
               "canonical_composition_L1\n";
        for (double pressureMPa : {26.0, 26.25, 26.5})
        {
            const double pressure = pressureMPa * 1.0e6;
            const Composition zHigh = compositionAt(parameters, 0.775);
            Flash::Result seed = flash.flashRestricted(
                pressure, temperature, zHigh, oilWaterPresence());
            if (!seed.converged)
                continue;
            Composition zTarget = failureZ;
            Flash::Result normal = flash.flashRestricted(
                pressure, temperature, zTarget,
                oilWaterPresence(), seed.composition);
            auto swappedSeed = seed.composition;
            std::swap(swappedSeed[0], swappedSeed[2]);
            Flash::Result swapped = flash.flashRestricted(
                pressure, temperature, zTarget,
                oilWaterPresence(), swappedSeed);
            role << pressureMPa << ','
                << (normal.converged ? 1 : 0) << ','
                << (swapped.converged ? 1 : 0) << ','
                << (normal.converged ? normal.composition[0][water] : missing) << ','
                << (normal.converged ? normal.composition[2][water] : missing) << ','
                << (swapped.converged ? swapped.composition[0][water] : missing) << ','
                << (swapped.converged ? swapped.composition[2][water] : missing) << ','
                << (normal.converged && swapped.converged
                    ? compositionDistance(normal, swapped) : missing)
                << '\n';
        }

        // Narrow direct grid to localize the failure manifold.
        std::ofstream grid(output / "narrow_direct_grid.csv");
        grid << std::setprecision(17)
            << "z_H2O,pressure_MPa,direct_converged,direct_phase_code,"
               "direct_iterations,water_only_stable,water_missing_oil_unstable,"
               "ow_unseeded_converged,ow_unseeded_stable\n";
        const std::array<double, 13> waterFractions{
            0.7400, 0.7450, 0.7500, 0.7550, 0.7600,
            failureWaterFraction, 0.7625, 0.7650, 0.7700,
            0.7750, 0.7800, 0.7850, 0.7900};
        for (double zw : waterFractions)
        {
            const Composition z = compositionAt(parameters, zw);
            for (int ip = 0; ip <= 60; ++ip)
            {
                const double pressureMPa = 25.5 + 0.025 * ip;
                const double pressure = pressureMPa * 1.0e6;
                const auto direct = certify(
                    eos, flash, pressure, z,
                    flash.flash(pressure, temperature, z));
                const auto waterOnly = certify(
                    eos, flash, pressure, z,
                    waterOnlyResult(flash, pressure, z));
                const auto ow = certify(
                    eos, flash, pressure, z,
                    flash.flashRestricted(
                        pressure, temperature, z, oilWaterPresence()));
                grid << zw << ',' << pressureMPa << ','
                    << (direct.converged ? 1 : 0) << ','
                    << direct.phaseCode << ',' << direct.iterations << ','
                    << (waterOnly.stabilityValid &&
                        waterOnly.stabilityStable ? 1 : 0) << ','
                    << (waterOnly.oilMissingUnstable ? 1 : 0) << ','
                    << (ow.converged ? 1 : 0) << ','
                    << (ow.stabilityValid && ow.stabilityStable ? 1 : 0)
                    << '\n';
            }
        }

        // Evidence-based classification at the two known failed states.
        std::ofstream classification(output / "classification.csv");
        classification
            << "hypothesis,status,evidence\n";

        bool seededBoth = true;
        bool unseededBoth = true;
        bool rootsBoth = true;
        bool roleBoth = true;
        bool stableBoth = true;
        for (double pressureMPa : {26.0, 26.25})
        {
            const double pressure = pressureMPa * 1.0e6;
            const auto unseeded = certify(
                eos, flash, pressure, failureZ,
                flash.flashRestricted(
                    pressure, temperature, failureZ, oilWaterPresence()));
            unseededBoth = unseededBoth && unseeded.converged;

            const Composition zHigh = compositionAt(parameters, 0.775);
            Flash::Result highSeed = flash.flashRestricted(
                pressure, temperature, zHigh, oilWaterPresence());
            Flash::Result seeded;
            if (highSeed.converged)
                seeded = flash.flashRestricted(
                    pressure, temperature, failureZ,
                    oilWaterPresence(), highSeed.composition);
            const auto seededCandidate =
                certify(eos, flash, pressure, failureZ, seeded);
            seededBoth = seededBoth && seededCandidate.converged;
            stableBoth = stableBoth &&
                seededCandidate.stabilityValid &&
                seededCandidate.stabilityStable;

            const auto probe = probeRoots(eos, pressure, failureZ);
            rootsBoth = rootsBoth &&
                probe.oilOk && probe.gasOk && probe.waterOk &&
                probe.profile.cpaAssociationIterativeCalls == 0 &&
                probe.profile.cpaWaterOnlyAnalyticFastPathAvailable;

            if (highSeed.converged)
            {
                auto swappedSeed = highSeed.composition;
                std::swap(swappedSeed[0], swappedSeed[2]);
                const auto normal = flash.flashRestricted(
                    pressure, temperature, failureZ,
                    oilWaterPresence(), highSeed.composition);
                const auto swapped = flash.flashRestricted(
                    pressure, temperature, failureZ,
                    oilWaterPresence(), swappedSeed);
                roleBoth = roleBoth && normal.converged && swapped.converged &&
                    compositionDistance(normal, swapped) < 1.0e-6;
            }
            else
                roleBoth = false;
        }

        classification
            << "CPA_ASSOCIATION_OR_DENSITY_ROOT_NUMERICS,"
            << (rootsBoth ? "NOT_SUPPORTED" : "POSSIBLE")
            << ",Direct Oil/Gas/Water CPA root probes plus association profiler at both failure states.\n"
            << "OIL_WATER_ROLE_CANONICALIZATION,"
            << (roleBoth ? "NOT_SUPPORTED" : "POSSIBLE")
            << ",Normal and deliberately swapped O/W seeds should converge to the same canonical tie-line.\n"
            << "UNSEEDED_ACTIVE_SET_OR_PAIR_INITIALIZATION,"
            << ((!unseededBoth && seededBoth) ? "SUPPORTED" : "NOT_ISOLATED")
            << ",Compare unseeded restricted O+W against the same solve with a nearby converged tie-line seed.\n"
            << "CONTINUATION_STEP_OR_SEED_REUSE,"
            << (seededBoth ? "SUPPORTED" : "NOT_RESOLVED")
            << ",A physically certified seeded O+W solution at both failed states shows a recoverable local equilibrium basin.\n"
            << "CPA_PARAMETER_OR_MODEL_FORM_PATHOLOGY,"
            << (seededBoth && stableBoth ? "NOT_REQUIRED_TO_EXPLAIN_FAILURE" : "POSSIBLE")
            << ",If a stable fugacity-closed O+W solution exists with frozen parameters, the direct failure is numerical rather than absence of an equilibrium solution.\n";

        std::cout << "CPA_ONSET_AUDIT_COMPLETE\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "CPA onset audit failed: " << error.what() << '\n';
        return 1;
    }
}
