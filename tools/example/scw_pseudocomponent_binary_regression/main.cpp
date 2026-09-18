/**
 * @file main.cpp
 * @brief Generic zero-flow H2O/pseudo-component PR regression and validation harness.
 *
 * This executable deliberately reuses the production CubicEquationOfState and
 * CubicThreePhaseFlash.  It fits only the PR binary interaction law to
 * coexistence-composition calibration rows.  Density/volume and phase-boundary
 * rows are independent validation gates and are never added to the BIP
 * objective, so kij cannot hide a density-model error.
 */
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <tools/eos_parameter_regression.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Config = MPMC::CompositionalModelConfig<
    2, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Mixture = MPMC::CompositionalMixture<Indices>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, 2>;
using Matrix = std::array<std::array<double, 2>, 2>;

constexpr double waterTc = 647.096;
constexpr double waterPc = 22.064e6;
constexpr double waterVc = 5.5948074534e-5;
constexpr double waterOmega = 0.3443;
constexpr double waterMw = 0.01801528;
constexpr double missing = std::numeric_limits<double>::quiet_NaN();

struct SystemDefinition
{
    std::string systemId;
    std::string proxyName;
    double hydrocarbonTc{};
    double hydrocarbonPc{};
    double hydrocarbonVc{};
    double hydrocarbonOmega{};
    double hydrocarbonMw{};
    double referenceTemperature{};
    double initialK{};
    double initialSlope{};
    double lowerK{};
    double upperK{};
    bool fitEnabled{false};
};

struct Observation
{
    std::string systemId;
    std::string id;
    std::string split;
    std::string observableType;
    double temperature{};
    double pressure{};
    double feedXWater{};
    int targetPhaseCount{};
    double targetXWaterRich{missing};
    double targetXHydrocarbonRich{missing};
    double toleranceXWaterRich{missing};
    double toleranceXHydrocarbonRich{missing};
    double targetDensityWaterRich{missing};
    double targetDensityHydrocarbonRich{missing};
    double densityRelativeTolerance{missing};
    std::string sourceId;
};

struct Prediction
{
    bool converged{false};
    int phaseCount{0};
    double xWaterRich{missing};
    double xHydrocarbonRich{missing};
    double densityWaterRich{missing};
    double densityHydrocarbonRich{missing};
};

std::vector<std::string> splitCsv(const std::string &line)
{
    std::vector<std::string> values;
    std::string current;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];
        if (c == '"')
        {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"')
            {
                current.push_back('"');
                ++i;
            }
            else
            {
                quoted = !quoted;
            }
        }
        else if (c == ',' && !quoted)
        {
            values.push_back(current);
            current.clear();
        }
        else
        {
            current.push_back(c);
        }
    }
    if (!current.empty() && current.back() == '\r')
        current.pop_back();
    values.push_back(current);
    return values;
}

std::map<std::string, std::size_t> readHeader(std::ifstream &input)
{
    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("CSV has no header");
    const auto fields = splitCsv(line);
    std::map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < fields.size(); ++i)
        columns.emplace(fields[i], i);
    return columns;
}

double parseNumber(const std::string &text)
{
    if (text.empty() || text == "nan" || text == "NA")
        return missing;
    return std::stod(text);
}

std::vector<SystemDefinition> readSystems(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open system manifest: " + path.string());
    const auto column = readHeader(input);
    std::vector<SystemDefinition> rows;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto value = splitCsv(line);
        auto text = [&](const std::string &name) -> const std::string & {
            return value.at(column.at(name));
        };
        rows.push_back({
            text("system_id"),
            text("proxy_name"),
            std::stod(text("Tc_K")),
            std::stod(text("Pc_MPa")) * 1.0e6,
            std::stod(text("Vc_cm3_mol")) * 1.0e-6,
            std::stod(text("omega")),
            std::stod(text("MW_g_mol")) * 1.0e-3,
            std::stod(text("T_ref_K")),
            std::stod(text("initial_kref")),
            std::stod(text("initial_b_K")),
            std::stod(text("k_lower")),
            std::stod(text("k_upper")),
            text("fit_enabled") == "1"});
    }
    if (rows.size() != 4)
        throw std::runtime_error("Expected exactly four pseudo-component systems");
    return rows;
}

std::vector<Observation> readObservations(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open observation data: " + path.string());
    const auto column = readHeader(input);
    std::vector<Observation> rows;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto value = splitCsv(line);
        auto text = [&](const std::string &name) -> const std::string & {
            return value.at(column.at(name));
        };
        rows.push_back({
            text("system_id"),
            text("id"),
            text("split"),
            text("observable_type"),
            std::stod(text("temperature_K")),
            std::stod(text("pressure_MPa")) * 1.0e6,
            std::stod(text("feed_x_water")),
            std::stoi(text("target_phase_count")),
            parseNumber(text("target_x_water_rich")),
            parseNumber(text("target_x_hydrocarbon_rich")),
            parseNumber(text("tolerance_x_water_rich")),
            parseNumber(text("tolerance_x_hydrocarbon_rich")),
            parseNumber(text("target_density_water_rich_kg_m3")),
            parseNumber(text("target_density_hydrocarbon_rich_kg_m3")),
            parseNumber(text("density_relative_tolerance")),
            text("source_id")});
    }
    return rows;
}

Mixture makeMixture(const SystemDefinition &system, double constantKij)
{
    const std::array<double, 2> tc{waterTc, system.hydrocarbonTc};
    const std::array<double, 2> pc{waterPc, system.hydrocarbonPc};
    const std::array<double, 2> vc{waterVc, system.hydrocarbonVc};
    const std::array<double, 2> omega{waterOmega, system.hydrocarbonOmega};
    const std::array<double, 2> mw{waterMw, system.hydrocarbonMw};
    const Matrix bip{{{{0.0, constantKij}}, {{constantKij, 0.0}}}};
    return Mixture(tc, pc, vc, omega, mw, bip);
}

Eos makeEos(
    const SystemDefinition &system,
    double kReference,
    double inverseTemperatureSlope)
{
    Eos eos(
        0.45724, 0.07780,
        makeMixture(system, kReference),
        1, 1.0 + std::sqrt(2.0), 1.0 - std::sqrt(2.0), 1.0e-30);
    eos.configureBinaryInteractionFunction(
        [system, kReference, inverseTemperatureSlope](
            int i, int j, double temperature) {
            if (i == j)
                return 0.0;
            return kReference + inverseTemperatureSlope *
                (1.0 / temperature - 1.0 / system.referenceTemperature);
        });
    return eos;
}

Prediction predict(
    const SystemDefinition &system,
    const Observation &row,
    double kReference,
    double inverseTemperatureSlope)
{
    const auto eos = makeEos(system, kReference, inverseTemperatureSlope);
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    const Flash flash(eos, options);
    const auto result = flash.flash(
        row.pressure, row.temperature,
        Composition{row.feedXWater, 1.0 - row.feedXWater});

    Prediction prediction;
    prediction.converged = result.converged;
    std::vector<std::size_t> active;
    if (result.converged)
    {
        for (std::size_t phase = 0; phase < 3; ++phase)
        {
            if (result.phaseMoleFraction[phase] > options.phaseFractionTolerance)
                active.push_back(phase);
        }
    }
    prediction.phaseCount = static_cast<int>(active.size());
    if (active.size() != 2)
        return prediction;

    std::sort(active.begin(), active.end(), [&](std::size_t a, std::size_t b) {
        return result.composition[a][0] > result.composition[b][0];
    });
    prediction.xWaterRich = result.composition[active[0]][0];
    prediction.xHydrocarbonRich = result.composition[active[1]][0];

    const auto density = [&](std::size_t slot) {
        const auto role = static_cast<MPMC::CompositionalPhase>(slot);
        const double molarDensity = eos.molarDensity(
            row.pressure, row.temperature, result.composition[slot],
            result.compressibility[slot], role);
        const double mixtureMw =
            result.composition[slot][0] * waterMw
            + result.composition[slot][1] * system.hydrocarbonMw;
        return molarDensity * mixtureMw;
    };
    prediction.densityWaterRich = density(active[0]);
    prediction.densityHydrocarbonRich = density(active[1]);
    return prediction;
}

bool hasCompositionTargets(const Observation &row)
{
    return std::isfinite(row.targetXWaterRich)
        && std::isfinite(row.targetXHydrocarbonRich);
}

double normalizedResidual(double predicted, double target, double tolerance)
{
    if (!(std::isfinite(tolerance) && tolerance > 0.0))
        throw std::runtime_error(
            "Experimental composition tolerance must be explicit and positive.");
    return (predicted - target) / tolerance;
}

double compositionObjective(
    const SystemDefinition &system,
    const std::vector<Observation> &rows,
    double kReference,
    double inverseTemperatureSlope)
{
    double sum = 0.0;
    int count = 0;
    try
    {
        for (const auto &row : rows)
        {
            if (row.systemId != system.systemId || row.split != "calibration"
                || row.observableType != "coexistence" || !hasCompositionTargets(row))
                continue;
            const double kij = kReference + inverseTemperatureSlope *
                (1.0 / row.temperature - 1.0 / system.referenceTemperature);
            if (kij < system.lowerK || kij > system.upperK)
                return 1.0e12;
            const auto prediction = predict(
                system, row, kReference, inverseTemperatureSlope);
            if (!prediction.converged || prediction.phaseCount != 2)
                return 1.0e12;
            const double rWater = normalizedResidual(
                prediction.xWaterRich, row.targetXWaterRich,
                row.toleranceXWaterRich);
            const double rHydrocarbon = normalizedResidual(
                prediction.xHydrocarbonRich, row.targetXHydrocarbonRich,
                row.toleranceXHydrocarbonRich);
            sum += 0.5 * (rWater * rWater + rHydrocarbon * rHydrocarbon);
            ++count;
        }
    }
    catch (const std::exception &)
    {
        return 1.0e12;
    }
    return count > 0 ? sum / count : 1.0e12;
}

struct FitState
{
    bool attempted{false};
    bool temperatureLaw{false};
    bool optimizerConverged{false};
    double kReference{};
    double slope{};
    double initialObjective{missing};
    double finalObjective{missing};
    int evaluations{0};
    std::string status;
};

FitState fitSystem(
    const SystemDefinition &system,
    const std::vector<Observation> &rows)
{
    FitState fit;
    fit.kReference = system.initialK;
    fit.slope = system.initialSlope;
    if (!system.fitEnabled)
    {
        fit.status = "DATA_OR_POLICY_BLOCKED";
        return fit;
    }

    std::set<double> calibrationTemperatures;
    int calibrationRows = 0;
    for (const auto &row : rows)
    {
        if (row.systemId == system.systemId && row.split == "calibration"
            && row.observableType == "coexistence" && hasCompositionTargets(row))
        {
            calibrationTemperatures.insert(row.temperature);
            ++calibrationRows;
        }
    }
    if (calibrationRows == 0)
    {
        fit.status = "BLOCKED_MISSING_RAW_ROWS";
        return fit;
    }

    fit.attempted = true;
    fit.temperatureLaw = calibrationTemperatures.size() >= 2;
    MPMC::tools::PatternSearchOptions options;
    options.maximumIterations = 180;
    options.minimumRelativeStep = 1.0e-9;
    options.improvementTolerance = 1.0e-14;

    std::vector<MPMC::tools::RegressionParameter> parameters{
        {"k_ref", system.initialK, system.lowerK, system.upperK, 0.01}};
    if (fit.temperatureLaw)
    {
        parameters.push_back(
            {"b_inverse_temperature_K", system.initialSlope,
             -5000.0, 5000.0, 25.0});
    }

    const auto result = MPMC::tools::boundedPatternSearch(
        parameters,
        [&](const std::vector<double> &x) {
            const double slope = fit.temperatureLaw ? x.at(1) : 0.0;
            return compositionObjective(system, rows, x.at(0), slope);
        }, options);

    fit.optimizerConverged = result.converged;
    fit.kReference = result.values.at(0);
    fit.slope = fit.temperatureLaw ? result.values.at(1) : 0.0;
    fit.initialObjective = result.initialObjective;
    fit.finalObjective = result.objective;
    fit.evaluations = result.objectiveEvaluations;
    fit.status = result.converged ? "REGRESSED_PENDING_PHYSICAL_GATES"
                                  : "REGRESSION_NOT_CONVERGED";
    return fit;
}

std::string passText(bool hasRows, bool passed)
{
    if (!hasRows)
        return "NO_DATA";
    return passed ? "PASS" : "FAIL";
}

struct GateState
{
    int phaseRows{0};
    int compositionRows{0};
    int densityRows{0};
    int boundaryRows{0};
    int validationRows{0};
    bool phasePass{true};
    bool compositionPass{true};
    bool densityPass{true};
    bool boundaryPass{true};
    bool holdoutPass{true};
};

bool densityWithin(double predicted, double target, double relativeTolerance)
{
    if (!std::isfinite(target))
        return true;
    if (!std::isfinite(predicted) || target <= 0.0)
        return false;
    if (!(std::isfinite(relativeTolerance) && relativeTolerance > 0.0))
        return false;
    return std::abs(predicted / target - 1.0) <= relativeTolerance;
}

GateState evaluateGates(
    const SystemDefinition &system,
    const std::vector<Observation> &rows,
    const FitState &fit,
    std::ofstream &predictions)
{
    GateState gate;
    for (const auto &row : rows)
    {
        if (row.systemId != system.systemId)
            continue;
        Prediction prediction;
        try
        {
            prediction = predict(system, row, fit.kReference, fit.slope);
        }
        catch (const std::exception &)
        {
            prediction = Prediction{};
        }

        bool phaseOk = true;
        if (row.targetPhaseCount > 0)
        {
            ++gate.phaseRows;
            phaseOk = prediction.converged
                && prediction.phaseCount == row.targetPhaseCount;
            gate.phasePass = gate.phasePass && phaseOk;
        }

        bool compositionOk = true;
        if (hasCompositionTargets(row))
        {
            ++gate.compositionRows;
            compositionOk = prediction.phaseCount == 2
                && std::abs(normalizedResidual(
                       prediction.xWaterRich, row.targetXWaterRich,
                       row.toleranceXWaterRich)) <= 1.0
                && std::abs(normalizedResidual(
                       prediction.xHydrocarbonRich,
                       row.targetXHydrocarbonRich,
                       row.toleranceXHydrocarbonRich)) <= 1.0;
            gate.compositionPass = gate.compositionPass && compositionOk;
        }

        const bool hasDensity = std::isfinite(row.targetDensityWaterRich)
            || std::isfinite(row.targetDensityHydrocarbonRich);
        bool densityOk = true;
        if (hasDensity)
        {
            ++gate.densityRows;
            densityOk = densityWithin(
                    prediction.densityWaterRich,
                    row.targetDensityWaterRich,
                    row.densityRelativeTolerance)
                && densityWithin(
                    prediction.densityHydrocarbonRich,
                    row.targetDensityHydrocarbonRich,
                    row.densityRelativeTolerance);
            gate.densityPass = gate.densityPass && densityOk;
        }

        bool boundaryOk = true;
        if (row.observableType == "boundary")
        {
            ++gate.boundaryRows;
            boundaryOk = phaseOk;
            gate.boundaryPass = gate.boundaryPass && boundaryOk;
        }

        if (row.split == "validation")
        {
            ++gate.validationRows;
            gate.holdoutPass = gate.holdoutPass
                && phaseOk && compositionOk && densityOk && boundaryOk;
        }

        const double kij = fit.kReference + fit.slope *
            (1.0 / row.temperature - 1.0 / system.referenceTemperature);
        predictions << system.systemId << ',' << system.proxyName << ','
            << row.id << ',' << row.split << ',' << row.observableType << ','
            << row.temperature << ',' << row.pressure / 1.0e6 << ',' << kij << ','
            << row.targetPhaseCount << ',' << prediction.phaseCount << ','
            << (prediction.converged ? 1 : 0) << ','
            << row.targetXWaterRich << ',' << prediction.xWaterRich << ','
            << row.targetXHydrocarbonRich << ',' << prediction.xHydrocarbonRich << ','
            << row.targetDensityWaterRich << ',' << prediction.densityWaterRich << ','
            << row.targetDensityHydrocarbonRich << ','
            << prediction.densityHydrocarbonRich << ',' << row.sourceId << '\n';
    }
    return gate;
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 4)
        {
            std::cerr
                << "Usage: scw_pseudocomponent_binary_regression "
                   "SYSTEMS.csv OBSERVATIONS.csv OUTPUT_DIR\n";
            return 2;
        }
        const auto systems = readSystems(argv[1]);
        const auto observations = readObservations(argv[2]);
        const std::filesystem::path output(argv[3]);
        if (std::filesystem::exists(output))
            throw std::runtime_error(
                "Output directory already exists: " + output.string());
        std::filesystem::create_directories(output);

        std::ofstream fits(output / "fit_parameters.csv");
        std::ofstream predictions(output / "predictions.csv");
        std::ofstream gates(output / "gate_summary.csv");
        if (!fits || !predictions || !gates)
            throw std::runtime_error("Cannot create regression output files");

        fits << std::setprecision(17)
            << "system_id,proxy_name,fit_status,temperature_law,optimizer_converged,"
               "k_ref,T_ref_K,b_inverse_temperature_K,initial_objective,final_objective,evaluations\n";
        predictions << std::setprecision(17)
            << "system_id,proxy_name,id,split,observable_type,temperature_K,pressure_MPa,kij,"
               "target_phase_count,predicted_phase_count,converged,target_x_water_rich,"
               "predicted_x_water_rich,target_x_hydrocarbon_rich,predicted_x_hydrocarbon_rich,"
               "target_density_water_rich_kg_m3,predicted_density_water_rich_kg_m3,"
               "target_density_hydrocarbon_rich_kg_m3,predicted_density_hydrocarbon_rich_kg_m3,source_id\n";
        gates << "system_id,fit_status,phase_rows,phase_gate,composition_rows,composition_gate,"
                 "density_rows,density_gate,boundary_rows,boundary_gate,validation_rows,holdout_gate,overall_status\n";

        bool allAccepted = true;
        for (const auto &system : systems)
        {
            const auto fit = fitSystem(system, observations);
            fits << system.systemId << ',' << system.proxyName << ',' << fit.status << ','
                << (fit.temperatureLaw ? 1 : 0) << ','
                << (fit.optimizerConverged ? 1 : 0) << ','
                << fit.kReference << ',' << system.referenceTemperature << ','
                << fit.slope << ',' << fit.initialObjective << ','
                << fit.finalObjective << ',' << fit.evaluations << '\n';

            const auto gate = evaluateGates(
                system, observations, fit, predictions);
            const bool accepted = fit.attempted && fit.optimizerConverged
                && gate.phaseRows > 0 && gate.phasePass
                && gate.compositionRows > 0 && gate.compositionPass
                && gate.densityRows > 0 && gate.densityPass
                && gate.boundaryRows > 0 && gate.boundaryPass
                && gate.validationRows > 0 && gate.holdoutPass;
            allAccepted = allAccepted && accepted;
            gates << system.systemId << ',' << fit.status << ','
                << gate.phaseRows << ',' << passText(gate.phaseRows > 0, gate.phasePass) << ','
                << gate.compositionRows << ','
                << passText(gate.compositionRows > 0, gate.compositionPass) << ','
                << gate.densityRows << ','
                << passText(gate.densityRows > 0, gate.densityPass) << ','
                << gate.boundaryRows << ','
                << passText(gate.boundaryRows > 0, gate.boundaryPass) << ','
                << gate.validationRows << ','
                << passText(gate.validationRows > 0, gate.holdoutPass) << ','
                << (accepted ? "PASS" : "BLOCKED") << '\n';
        }

        std::ofstream overall(output / "reservoir_gate.txt");
        overall << (allAccepted ? "RESERVOIR_GATE_PASS\n"
                               : "RESERVOIR_GATE_BLOCKED\n");
        std::cout << (allAccepted ? "RESERVOIR_GATE_PASS\n"
                                 : "RESERVOIR_GATE_BLOCKED\n");
        return allAccepted ? 0 : 3;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
