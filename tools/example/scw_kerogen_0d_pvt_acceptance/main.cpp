/**
 * @file main.cpp
 * @brief Zero-flow PVT acceptance for the experiment-driven SCW kerogen fluid.
 *
 * PR and CPA are evaluated independently on the same registered T/P/z grid.
 * This executable reuses production CubicEquationOfState,
 * CubicThreePhaseFlash, stabilityTest and phase-diagram utilities.
 */
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/properties/aqueous_viscosity.hpp>
#include <natural/properties/compositional_properties.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <tools/phase_diagram.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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
using PropertyModel = MPMC::CompositionalPropertyModel<Indices>;
using Composition = std::array<double, 5>;
using Matrix = std::array<std::array<double, 5>, 5>;

constexpr std::size_t N = 5;
constexpr std::size_t water = 0;
constexpr double tref = 653.15;
constexpr double missing = std::numeric_limits<double>::quiet_NaN();

inline const std::array<std::string, N> names{
    "H2O", "OIL_GASOLINE", "OIL_DIESEL", "OIL_MIDDLE", "OIL_HEAVY"};
inline const std::array<const char *, 3> phaseRoleNames{"Oil", "Gas", "Water"};
inline constexpr std::array<double, 3> targetTemperatures{
    633.15, 647.15, 653.15};
inline constexpr std::array<double, 6> targetPressuresMPa{
    25.0, 26.0, 27.0, 28.0, 29.0, 30.0};

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
    Matrix prKref{};
    Matrix prSlope{};
    std::array<double, N> cpaA0{};
    std::array<double, N> cpaB{};
    std::array<double, N> cpaC1{};
    std::array<double, N> cpaAssociationEnergy{};
    std::array<double, N> cpaAssociationVolume{};
    std::array<int, N> cpaDonorSites{};
    std::array<int, N> cpaAcceptorSites{};
    Matrix cpaKij{};
};

struct ScanComposition
{
    std::string family;
    double waterFraction{};
    bool initialAnchor{false};
    Composition z{};
};

struct PhaseProperty
{
    bool active{false};
    double beta{missing};
    double saturation{missing};
    double zFactor{missing};
    double molarDensity{missing};
    double massDensity{missing};
    double lbcViscosity{missing};
    double iapwsViscosity{missing};
};

struct StateEvaluation
{
    bool flashConverged{false};
    int phaseCode{0};
    int phaseCount{0};
    bool stabilityValid{false};
    bool stabilityStable{false};
    std::array<double, 3> trialSum{missing, missing, missing};
    std::array<bool, 3> missingPhaseUnstable{false, false, false};
    double massClosure{missing};
    double maxLogFugacitySpread{missing};
    bool waterRoleConsistent{false};
    bool phaseCompositionsValid{false};
    bool phasePropertiesFinite{false};
    bool pass{false};
    Flash::Result flash{};
    std::array<PhaseProperty, 3> property{};
};

struct AnchorResult
{
    double temperature{};
    double pressureMPa{};
    int phaseCode{};
    int phaseCount{};
    bool pass{false};
};

struct BackendResult
{
    std::string name;
    std::size_t scanPoints{0};
    std::size_t passedPoints{0};
    std::size_t initialPoints{0};
    std::size_t passedInitialPoints{0};
    std::vector<AnchorResult> anchors;
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
    {
        if (!line.empty())
            result.row.push_back(splitCsv(line));
    }
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
    throw std::runtime_error("Unknown component name: " + name);
}

void normalize(Composition &z)
{
    double sum = 0.0;
    for (double value : z)
    {
        if (!std::isfinite(value) || value < 0.0)
            throw std::runtime_error("Composition is not finite/non-negative.");
        sum += value;
    }
    if (!(sum > 0.0))
        throw std::runtime_error("Composition sum is not positive.");
    for (double &value : z)
        value /= sum;
}

Parameters readParameters(const std::filesystem::path &caseDir)
{
    Parameters p;

    const auto prPure = readCsv(
        caseDir / "pr_parameters/pr_pure_parameters_380c.csv");
    std::size_t pureRows = 0;
    for (const auto &row : prPure.row)
    {
        const std::size_t i = componentIndex(field(prPure, row, "component"));
        p.tc[i] = std::stod(field(prPure, row, "Tc_K"));
        p.pc[i] = std::stod(field(prPure, row, "Pc_MPa")) * 1.0e6;
        p.vc[i] = std::stod(field(prPure, row, "Vc_cm3_mol")) * 1.0e-6;
        p.omega[i] = std::stod(field(prPure, row, "omega"));
        p.mw[i] = std::stod(field(prPure, row, "MW_g_mol")) * 1.0e-3;
        ++pureRows;
    }
    if (pureRows != N)
        throw std::runtime_error("Expected five PR pure-component rows.");

    const auto prBip = readCsv(
        caseDir / "pr_parameters/pr_binary_matrix_screening.csv");
    for (const auto &row : prBip.row)
    {
        const std::size_t i = componentIndex(field(prBip, row, "component_i"));
        const std::size_t j = componentIndex(field(prBip, row, "component_j"));
        const double k = std::stod(field(prBip, row, "kref_at_65315K"));
        const double b = std::stod(field(prBip, row, "b_inverse_temperature_K"));
        p.prKref[i][j] = p.prKref[j][i] = k;
        p.prSlope[i][j] = p.prSlope[j][i] = b;
    }

    const auto cpaPure = readCsv(
        caseDir / "cpa_parameters/cpa_pure_parameters_380c.csv");
    std::size_t cpaRows = 0;
    for (const auto &row : cpaPure.row)
    {
        const std::size_t i = componentIndex(field(cpaPure, row, "component"));
        p.cpaA0[i] = std::stod(field(cpaPure, row, "cpa_a0_Pa_m6_mol2"));
        p.cpaB[i] = std::stod(field(cpaPure, row, "cpa_b_m3_mol"));
        p.cpaC1[i] = std::stod(field(cpaPure, row, "cpa_c1"));
        ++cpaRows;
    }
    if (cpaRows != N)
        throw std::runtime_error("Expected five CPA pure-component rows.");

    const auto assoc = readCsv(
        caseDir / "cpa_parameters/cpa_association_scheme.csv");
    for (const auto &row : assoc.row)
    {
        const std::size_t i = componentIndex(field(assoc, row, "component"));
        p.cpaDonorSites[i] = std::stoi(field(assoc, row, "donor_sites"));
        p.cpaAcceptorSites[i] = std::stoi(field(assoc, row, "acceptor_sites"));
        p.cpaAssociationEnergy[i] =
            std::stod(field(assoc, row, "association_energy_J_mol"));
        p.cpaAssociationVolume[i] =
            std::stod(field(assoc, row, "association_volume"));
    }

    const auto cpaBip = readCsv(
        caseDir / "cpa_parameters/cpa_binary_matrix_screening.csv");
    for (const auto &row : cpaBip.row)
    {
        const std::size_t i = componentIndex(field(cpaBip, row, "component_i"));
        const std::size_t j = componentIndex(field(cpaBip, row, "component_j"));
        const double k = std::stod(field(cpaBip, row, "kij_constant_initial"));
        p.cpaKij[i][j] = p.cpaKij[j][i] = k;
    }
    return p;
}

std::vector<ScanComposition> readScan(
    const std::filesystem::path &caseDir)
{
    const auto csv = readCsv(caseDir / "pvt_acceptance/composition_scan.csv");
    std::vector<ScanComposition> result;
    for (const auto &row : csv.row)
    {
        ScanComposition item;
        item.family = field(csv, row, "family");
        item.waterFraction = std::stod(field(csv, row, "water_fraction"));
        item.initialAnchor = field(csv, row, "initial_anchor") == "1";
        item.z = {
            std::stod(field(csv, row, "z_H2O")),
            std::stod(field(csv, row, "z_OIL_GASOLINE")),
            std::stod(field(csv, row, "z_OIL_DIESEL")),
            std::stod(field(csv, row, "z_OIL_MIDDLE")),
            std::stod(field(csv, row, "z_OIL_HEAVY"))};
        normalize(item.z);
        result.push_back(item);
    }
    return result;
}

Mixture makeMixture(
    const Parameters &p,
    const Matrix &kij)
{
    return Mixture(p.tc, p.pc, p.vc, p.omega, p.mw, kij);
}

Eos makePr(const Parameters &p)
{
    Eos eos(
        0.45724, 0.07780, makeMixture(p, p.prKref),
        1, 1.0 + std::sqrt(2.0), 1.0 - std::sqrt(2.0), 1.0e-30);
    const Matrix kref = p.prKref;
    const Matrix slope = p.prSlope;
    eos.configureBinaryInteractionFunction(
        [kref, slope](int i, int j, double temperature) {
            if (i == j)
                return 0.0;
            const auto ii = static_cast<std::size_t>(i);
            const auto jj = static_cast<std::size_t>(j);
            return kref[ii][jj] + slope[ii][jj]
                * (1.0 / temperature - 1.0 / tref);
        });
    return eos;
}

Eos makeCpa(const Parameters &p)
{
    Eos eos(
        0.42748, 0.08664, makeMixture(p, p.cpaKij),
        1, 1.0, 0.0, 1.0e-30);
    Eos::CubicPlusAssociationOptions cpa;
    cpa.a0 = p.cpaA0;
    cpa.b = p.cpaB;
    cpa.c1 = p.cpaC1;
    cpa.associationEnergy = p.cpaAssociationEnergy;
    cpa.associationVolume = p.cpaAssociationVolume;
    cpa.donorSites = p.cpaDonorSites;
    cpa.acceptorSites = p.cpaAcceptorSites;
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::SoaveRedlichKwong;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

double massClosure(const Composition &z, const Flash::Result &result)
{
    double maximum = 0.0;
    for (std::size_t component = 0; component < N; ++component)
    {
        double reconstructed = 0.0;
        for (std::size_t phase = 0; phase < 3; ++phase)
            reconstructed += result.phaseMoleFraction[phase]
                * result.composition[phase][component];
        maximum = std::max(
            maximum, std::abs(reconstructed - z[component]));
    }
    return maximum;
}

double maxLogFugacitySpread(
    const Eos &eos,
    double pressure,
    double temperature,
    const Flash::Result &result)
{
    if (result.presence.count() <= 1)
        return 0.0;
    std::array<Eos::PhaseResult<double>, 3> thermo{};
    std::array<bool, 3> active{false, false, false};
    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        const auto role = static_cast<MPMC::CompositionalPhase>(phase);
        if (!result.presence.contains(role))
            continue;
        active[phase] = true;
        thermo[phase] = eos.phaseResult(
            pressure, temperature, result.composition[phase], role, false);
    }

    double maximum = 0.0;
    for (std::size_t component = 0; component < N; ++component)
    {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        int count = 0;
        for (std::size_t phase = 0; phase < 3; ++phase)
        {
            if (!active[phase])
                continue;
            const double fugacity = std::max(
                thermo[phase].fugacity[component], 1.0e-300);
            const double logF = std::log(fugacity);
            lo = std::min(lo, logF);
            hi = std::max(hi, logF);
            ++count;
        }
        if (count > 1)
            maximum = std::max(maximum, hi - lo);
    }
    return maximum;
}

bool phaseCompositionsValid(const Flash::Result &result)
{
    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        const auto role = static_cast<MPMC::CompositionalPhase>(phase);
        if (!result.presence.contains(role))
            continue;
        double sum = 0.0;
        for (double value : result.composition[phase])
        {
            if (!std::isfinite(value) || value < -1.0e-12)
                return false;
            sum += value;
        }
        if (std::abs(sum - 1.0) > 1.0e-8)
            return false;
    }
    return true;
}

bool waterRoleConsistent(const Flash::Result &result)
{
    if (!result.presence.contains(MPMC::CompositionalPhase::Water))
        return true;
    const double waterRole =
        result.composition[MPMC::phaseIndex(MPMC::CompositionalPhase::Water)][water];
    for (std::size_t phase = 0; phase < 2; ++phase)
    {
        const auto role = static_cast<MPMC::CompositionalPhase>(phase);
        if (result.presence.contains(role) &&
            waterRole + 1.0e-8 < result.composition[phase][water])
            return false;
    }
    return true;
}

StateEvaluation evaluateState(
    const Eos &eos,
    const Flash &flash,
    double pressure,
    double temperature,
    Composition z)
{
    normalize(z);
    StateEvaluation state;
    state.flash = flash.flash(pressure, temperature, z);
    state.flashConverged = state.flash.converged;
    if (!state.flashConverged)
        return state;

    state.phaseCode = static_cast<int>(state.flash.presence.bits());
    state.phaseCount = state.flash.presence.count();
    const auto stability = flash.stabilityTest(
        pressure, temperature, z,
        state.flash.presence, state.flash.composition);
    state.stabilityValid = stability.valid;
    state.stabilityStable = stability.stable;
    state.trialSum = stability.trialSum;
    state.missingPhaseUnstable = stability.missingPhaseUnstable;
    state.massClosure = massClosure(z, state.flash);
    state.maxLogFugacitySpread =
        maxLogFugacitySpread(eos, pressure, temperature, state.flash);
    state.waterRoleConsistent = waterRoleConsistent(state.flash);
    state.phaseCompositionsValid = phaseCompositionsValid(state.flash);

    const PropertyModel properties(eos.mixture());
    const MPMC::Iapws2008IndustrialAqueousViscosity<5> iapws(
        static_cast<int>(water), 0.02);
    state.phasePropertiesFinite = true;
    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        auto &property = state.property[phase];
        const auto role = static_cast<MPMC::CompositionalPhase>(phase);
        property.active = state.flash.presence.contains(role);
        if (!property.active)
            continue;
        property.beta = state.flash.phaseMoleFraction[phase];
        property.saturation = state.flash.saturation[phase];
        property.zFactor = state.flash.compressibility[phase];
        property.molarDensity = state.flash.molarDensity[phase];
        double phaseMw = 0.0;
        for (std::size_t component = 0; component < N; ++component)
            phaseMw += state.flash.composition[phase][component]
                * eos.mixture().molecularWeight(static_cast<int>(component));
        property.massDensity = property.molarDensity * phaseMw;
        property.lbcViscosity = properties.viscosityFromMolarDensity(
            state.flash.composition[phase], property.molarDensity, temperature);
        if (role == MPMC::CompositionalPhase::Water)
        {
            try
            {
                property.iapwsViscosity = iapws.viscosity(
                    pressure, temperature, state.flash.composition[phase]);
            }
            catch (const std::exception &)
            {
                property.iapwsViscosity = missing;
            }
        }

        const bool finite =
            std::isfinite(property.beta) && property.beta >= 0.0 &&
            std::isfinite(property.saturation) && property.saturation >= 0.0 &&
            std::isfinite(property.zFactor) && property.zFactor > 0.0 &&
            std::isfinite(property.molarDensity) && property.molarDensity > 0.0 &&
            std::isfinite(property.massDensity) && property.massDensity > 0.0 &&
            std::isfinite(property.lbcViscosity) && property.lbcViscosity > 0.0;
        state.phasePropertiesFinite = state.phasePropertiesFinite && finite;
    }

    state.pass =
        state.phaseCount >= 1 &&
        state.stabilityValid &&
        state.stabilityStable &&
        state.massClosure <= 1.0e-8 &&
        state.maxLogFugacitySpread <= 1.0e-6 &&
        state.waterRoleConsistent &&
        state.phaseCompositionsValid &&
        state.phasePropertiesFinite;
    return state;
}

std::string safeStem(std::string name)
{
    for (char &ch : name)
        if (!(std::isalnum(static_cast<unsigned char>(ch)) || ch == '_'))
            ch = '_';
    return name;
}

BackendResult runBackend(
    const std::string &backend,
    Eos eos,
    const std::vector<ScanComposition> &scan,
    const std::filesystem::path &outputRoot)
{
    const std::filesystem::path outDir = outputRoot / backend;
    std::filesystem::create_directories(outDir);

    MPMC::ThreePhaseFlashOptions flashOptions;
    flashOptions.waterComponent = static_cast<int>(water);
    flashOptions.maximumIterations = 180;
    flashOptions.maximumStabilityIterations = 120;
    const Flash flash(eos, flashOptions);

    std::ofstream states(outDir / "state_scan.csv");
    std::ofstream phases(outDir / "phase_properties.csv");
    if (!states || !phases)
        throw std::runtime_error("Cannot create 0D PVT output files.");

    states << std::setprecision(17)
        << "backend,family,water_fraction,initial_anchor,temperature_C,temperature_K,"
           "pressure_MPa,phase_code,phase_count,converged,stability_valid,stability_stable,"
           "trial_sum_o,trial_sum_g,trial_sum_w,missing_unstable_o,missing_unstable_g,"
           "missing_unstable_w,mass_closure,max_log_fugacity_spread,water_role_consistent,"
           "phase_compositions_valid,phase_properties_finite,state_pass\n";
    phases << std::setprecision(17)
        << "backend,family,water_fraction,initial_anchor,temperature_C,temperature_K,"
           "pressure_MPa,phase_role,active,beta,saturation,Z,molar_density_mol_m3,"
           "mass_density_kg_m3,viscosity_LBC_Pa_s,viscosity_IAPWS_if_admissible_Pa_s";
    for (const auto &name : names)
        phases << ",x_" << name;
    phases << '\n';

    BackendResult summary;
    summary.name = backend;
    for (const auto &composition : scan)
    {
        for (double temperature : targetTemperatures)
        {
            for (double pressureMPa : targetPressuresMPa)
            {
                const auto state = evaluateState(
                    eos, flash, pressureMPa * 1.0e6,
                    temperature, composition.z);
                ++summary.scanPoints;
                summary.passedPoints += state.pass ? 1u : 0u;
                const bool anchor = composition.initialAnchor
                    && std::abs(pressureMPa - 25.0) < 1.0e-12;
                if (anchor)
                {
                    ++summary.initialPoints;
                    summary.passedInitialPoints += state.pass ? 1u : 0u;
                    summary.anchors.push_back({
                        temperature, pressureMPa,
                        state.phaseCode, state.phaseCount, state.pass});
                }

                states << backend << ',' << composition.family << ','
                    << composition.waterFraction << ','
                    << (composition.initialAnchor ? 1 : 0) << ','
                    << temperature - 273.15 << ',' << temperature << ','
                    << pressureMPa << ',' << state.phaseCode << ','
                    << state.phaseCount << ','
                    << (state.flashConverged ? 1 : 0) << ','
                    << (state.stabilityValid ? 1 : 0) << ','
                    << (state.stabilityStable ? 1 : 0);
                for (double value : state.trialSum)
                    states << ',' << value;
                for (bool value : state.missingPhaseUnstable)
                    states << ',' << (value ? 1 : 0);
                states << ',' << state.massClosure
                    << ',' << state.maxLogFugacitySpread
                    << ',' << (state.waterRoleConsistent ? 1 : 0)
                    << ',' << (state.phaseCompositionsValid ? 1 : 0)
                    << ',' << (state.phasePropertiesFinite ? 1 : 0)
                    << ',' << (state.pass ? 1 : 0) << '\n';

                for (std::size_t phase = 0; phase < 3; ++phase)
                {
                    const auto &property = state.property[phase];
                    phases << backend << ',' << composition.family << ','
                        << composition.waterFraction << ','
                        << (composition.initialAnchor ? 1 : 0) << ','
                        << temperature - 273.15 << ',' << temperature << ','
                        << pressureMPa << ',' << phaseRoleNames[phase] << ','
                        << (property.active ? 1 : 0) << ','
                        << property.beta << ',' << property.saturation << ','
                        << property.zFactor << ',' << property.molarDensity << ','
                        << property.massDensity << ',' << property.lbcViscosity << ','
                        << property.iapwsViscosity;
                    for (double x : state.flash.composition[phase])
                        phases << ',' << x;
                    phases << '\n';
                }
            }
        }
    }

    // Full O/G/W target-window maps and phase-onset boundaries at z_H2O=0.20
    // for each deterministic oil-composition family.
    MPMC::tools::PhaseDiagramSampler<Indices> unrestricted(
        flash, MPMC::tools::PhaseDiagramFlashPolicy::unrestricted());
    const auto oilGas = MPMC::PhasePresence(
        static_cast<std::uint8_t>(
            MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit));
    MPMC::tools::PhaseDiagramSampler<Indices> restricted(
        flash, MPMC::tools::PhaseDiagramFlashPolicy::restricted(oilGas));
    const MPMC::tools::ScanAxis envelopeT{
        628.15, 658.15, 31, MPMC::tools::AxisSpacing::Linear};
    const MPMC::tools::ScanAxis envelopeP{
        20.0e6, 35.0e6, 61, MPMC::tools::AxisSpacing::Linear};
    MPMC::tools::MultiphaseBoundaryOptions boundaryOptions;
    boundaryOptions.maxRefinementIterations = 24;
    boundaryOptions.relativePressureTolerance = 1.0e-6;
    MPMC::tools::EnvelopeOptions envelopeOptions;
    envelopeOptions.maxRefinementIterations = 24;
    envelopeOptions.relativePressureTolerance = 1.0e-6;

    for (const std::string family : {"BASE", "LIGHT_ENRICHED", "HEAVY_ENRICHED"})
    {
        auto it = std::find_if(
            scan.begin(), scan.end(), [&](const ScanComposition &item) {
                return item.family == family &&
                    std::abs(item.waterFraction - 0.20) < 1.0e-12;
            });
        if (it == scan.end())
            throw std::runtime_error("Missing z_H2O=0.20 family scan row.");
        const std::string stem = safeStem(family);
        const auto pt = unrestricted.pressureTemperature(
            it->z, envelopeT, envelopeP);
        MPMC::tools::writeCsv(
            outDir / (stem + "_full_pt_map.csv"), pt, names);
        const auto boundaries =
            unrestricted.pressureTemperaturePhaseBoundaries(
                pt, boundaryOptions);
        MPMC::tools::writeCsv(
            outDir / (stem + "_phase_onset_envelope.csv"), boundaries);
        const auto ogEnvelope = restricted.pressureTemperatureEnvelope(
            it->z, envelopeT, envelopeP, envelopeOptions);
        MPMC::tools::writeCsv(
            outDir / (stem + "_restricted_og_envelope.csv"),
            ogEnvelope, names);
    }

    // Dense pressure-composition scan along the measured oil ratio for each
    // requested target temperature.
    auto low = std::find_if(
        scan.begin(), scan.end(), [](const ScanComposition &item) {
            return item.family == "BASE" &&
                std::abs(item.waterFraction - 0.01) < 1.0e-12;
        });
    auto high = std::find_if(
        scan.begin(), scan.end(), [](const ScanComposition &item) {
            return item.family == "BASE" &&
                std::abs(item.waterFraction - 0.995) < 1.0e-12;
        });
    if (low == scan.end() || high == scan.end())
        throw std::runtime_error("Missing BASE composition-path endpoints.");
    const MPMC::tools::ScanAxis targetPressure{
        25.0e6, 30.0e6, 21, MPMC::tools::AxisSpacing::Linear};
    const MPMC::tools::ScanAxis pathAxis{
        0.0, 1.0, 81, MPMC::tools::AxisSpacing::Linear};
    for (double temperature : targetTemperatures)
    {
        const auto map = unrestricted.pressureComposition(
            temperature, targetPressure, low->z, high->z, pathAxis);
        std::ostringstream stem;
        stem << "base_pressure_composition_" << std::llround(temperature * 100.0);
        MPMC::tools::writeCsv(outDir / (stem.str() + ".csv"), map, names);
    }

    std::ofstream summaryFile(outDir / "acceptance_summary.csv");
    summaryFile << "backend,scan_points,scan_passed,scan_failed,scan_pass_fraction,"
                   "initial_anchor_points,initial_anchor_passed,initial_anchor_failed,"
                   "scan_health,initial_state_gate\n"
        << backend << ',' << summary.scanPoints << ',' << summary.passedPoints << ','
        << (summary.scanPoints - summary.passedPoints) << ','
        << static_cast<double>(summary.passedPoints) /
               static_cast<double>(summary.scanPoints) << ','
        << summary.initialPoints << ',' << summary.passedInitialPoints << ','
        << (summary.initialPoints - summary.passedInitialPoints) << ','
        << (summary.passedPoints == summary.scanPoints ? "PASS" : "FAIL") << ','
        << (summary.passedInitialPoints == summary.initialPoints ? "PASS" : "FAIL")
        << '\n';
    return summary;
}

const AnchorResult *findAnchor(
    const BackendResult &result,
    double temperature,
    double pressureMPa)
{
    for (const auto &anchor : result.anchors)
        if (std::abs(anchor.temperature - temperature) < 1.0e-10 &&
            std::abs(anchor.pressureMPa - pressureMPa) < 1.0e-10)
            return &anchor;
    return nullptr;
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            std::cerr
                << "Usage: scw_kerogen_0d_pvt_acceptance CASE_DIR OUTPUT_DIR\n";
            return 2;
        }
        const std::filesystem::path caseDir(argv[1]);
        const std::filesystem::path output(argv[2]);
        if (std::filesystem::exists(output))
            throw std::runtime_error(
                "Output directory already exists: " + output.string());
        std::filesystem::create_directories(output);

        const Parameters parameters = readParameters(caseDir);
        const auto scan = readScan(caseDir);

        const auto pr = runBackend(
            "PR", makePr(parameters), scan, output);
        const auto cpa = runBackend(
            "CPA", makeCpa(parameters), scan, output);

        std::ofstream cross(output / "cross_eos_initial_state.csv");
        cross << "temperature_C,temperature_K,pressure_MPa,"
                 "pr_phase_code,pr_phase_count,pr_state_pass,"
                 "cpa_phase_code,cpa_phase_count,cpa_state_pass,"
                 "both_physically_self_consistent,phase_count_agreement_diagnostic\n";
        bool bothInitial = true;
        for (double temperature : targetTemperatures)
        {
            const auto *prAnchor = findAnchor(pr, temperature, 25.0);
            const auto *cpaAnchor = findAnchor(cpa, temperature, 25.0);
            if (!prAnchor || !cpaAnchor)
                throw std::runtime_error("Missing initial-state anchor result.");
            const bool both = prAnchor->pass && cpaAnchor->pass;
            bothInitial = bothInitial && both;
            cross << temperature - 273.15 << ',' << temperature << ",25,"
                << prAnchor->phaseCode << ',' << prAnchor->phaseCount << ','
                << (prAnchor->pass ? 1 : 0) << ','
                << cpaAnchor->phaseCode << ',' << cpaAnchor->phaseCount << ','
                << (cpaAnchor->pass ? 1 : 0) << ','
                << (both ? 1 : 0) << ','
                << (prAnchor->phaseCount == cpaAnchor->phaseCount ? 1 : 0)
                << '\n';
        }

        const bool prScan =
            pr.scanPoints > 0 && pr.passedPoints == pr.scanPoints;
        const bool cpaScan =
            cpa.scanPoints > 0 && cpa.passedPoints == cpa.scanPoints;
        const bool zeroDPvtPass = prScan && cpaScan && bothInitial;

        std::ofstream gate(output / "zero_d_pvt_gate.csv");
        gate << "gate,status,requirement\n"
            << "PR_TARGET_SCAN," << (prScan ? "PASS" : "FAIL")
            << ",All registered PR T-P-z scan states must be physically self-consistent\n"
            << "CPA_TARGET_SCAN," << (cpaScan ? "PASS" : "FAIL")
            << ",All registered CPA T-P-z scan states must be physically self-consistent\n"
            << "CROSS_EOS_INITIAL_STATE," << (bothInitial ? "PASS" : "FAIL")
            << ",At 360/374/380 C and 25 MPa the BASE initial composition must be self-consistent in both EOS; identical phase count is diagnostic only\n"
            << "ZERO_D_PVT_ACCEPTANCE," << (zeroDPvtPass ? "PASS" : "BLOCKED")
            << ",PR and CPA target scans plus cross-EOS initial-state gate\n";

        std::ofstream flag(output / "zero_d_pvt_gate.txt");
        flag << (zeroDPvtPass
            ? "ZERO_D_PVT_GATE_PASS\n"
            : "ZERO_D_PVT_GATE_BLOCKED\n");

        std::cout << (zeroDPvtPass
            ? "ZERO_D_PVT_GATE_PASS\n"
            : "ZERO_D_PVT_GATE_BLOCKED\n");
        return zeroDPvtPass ? 0 : 3;
    }
    catch (const std::exception &error)
    {
        std::cerr << "0D PVT acceptance failed: " << error.what() << '\n';
        return 1;
    }
}
