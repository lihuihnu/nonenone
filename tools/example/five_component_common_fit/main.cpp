/**
 * @file main.cpp
 * @brief 使用固定公开数据集拟合 PR、SW 与 CPA 的组分对参数。
 */
#include "../../../case/five_component_eos_tuned_compare/case_config.hpp"

#include <indices/indices.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/soreide_whitson.hpp>
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
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{

using FitConfig = MPMC::CompositionalModelConfig<
    5, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<FitConfig>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Composition = std::array<double, 5>;
using Matrix = std::array<std::array<double, 5>, 5>;

constexpr std::size_t h2o = 0;
constexpr std::size_t co2 = 1;
constexpr std::size_t ch4 = 2;
constexpr std::size_t c2 = 3;
constexpr std::size_t nc4 = 4;

struct PairDefinition
{
    std::string name;
    std::size_t first;
    std::size_t second;
    bool aqueous;
    double lower;
    double upper;
    double step;
};

const std::array<PairDefinition, 6> pairs{{
    {"h2o_co2", h2o, co2, true, -0.20, 0.60, 0.04},
    {"h2o_ch4", h2o, ch4, true, -0.20, 0.65, 0.04},
    {"co2_ch4", co2, ch4, false, -0.20, 0.30, 0.02},
    {"co2_c2", co2, c2, false, -0.20, 0.30, 0.02},
    {"ch4_c2", ch4, c2, false, -0.20, 0.30, 0.02},
    {"ch4_nc4", ch4, nc4, false, -0.20, 0.30, 0.02},
}};

struct Datum
{
    std::string id;
    std::string source;
    std::string pair;
    std::string kind;
    std::string split;
    double temperature{0.0};
    double pressure{0.0};
    double salinity{0.0};
    double density{0.0};
    double uncertainty{0.0};
    Composition liquid{};
    Composition vapor{};
};

struct PointPrediction
{
    std::string quantity;
    std::string component;
    double observed{0.0};
    double predicted{0.0};
    double relativeError{0.0};
};

struct FitRecord
{
    std::string model;
    std::string pair;
    std::string parameterKind;
    double initial{0.0};
    double lower{0.0};
    double upper{0.0};
    double fitted{0.0};
    double trainingObjectiveInitial{0.0};
    double trainingObjectiveFinal{0.0};
    int evaluations{0};
};

std::vector<std::string> splitCsv(const std::string &line)
{
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
    {
        if (!value.empty() && value.back() == '\r')
            value.pop_back();
        values.push_back(value);
    }
    if (!line.empty() && line.back() == ',')
        values.emplace_back();
    return values;
}

std::vector<Datum> readData(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open common calibration CSV: " + path.string());
    std::string line;
    std::getline(input, line);
    const auto header = splitCsv(line);
    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t index = 0; index < header.size(); ++index)
        column.emplace(header[index], index);
    auto text = [&](const std::vector<std::string> &row, const std::string &name) {
        return row.at(column.at(name));
    };
    auto number = [&](const std::vector<std::string> &row, const std::string &name) {
        return std::stod(text(row, name));
    };

    std::vector<Datum> data;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto row = splitCsv(line);
        Datum datum;
        datum.id = text(row, "id");
        datum.source = text(row, "source");
        datum.pair = text(row, "pair");
        datum.kind = text(row, "kind");
        datum.split = text(row, "split");
        datum.temperature = number(row, "temperature_K");
        datum.pressure = number(row, "pressure_Pa");
        datum.salinity = number(row, "salinity_molal");
        datum.density = number(row, "density_kg_m3");
        datum.uncertainty = number(row, "uncertainty");
        const std::array<std::string, 5> names{"h2o", "co2", "ch4", "c2", "nc4"};
        for (std::size_t component = 0; component < names.size(); ++component)
        {
            datum.liquid[component] = number(row, "x_" + names[component]);
            datum.vapor[component] = number(row, "y_" + names[component]);
        }
        data.push_back(std::move(datum));
    }
    return data;
}

Matrix baseMatrix(const std::string &model)
{
    if (model == "pr")
        return CaseConfig::PrFluid::binaryInteraction;
    if (model == "sw")
        return CaseConfig::SwFluid::binaryInteraction;
    if (model == "cpa")
        return CaseConfig::CpaFluid::binaryInteraction;
    throw std::invalid_argument("Unknown EOS model: " + model);
}

std::array<double, 6> initialParameters(const std::string &model)
{
    const auto matrix = baseMatrix(model);
    std::array<double, 6> values{};
    for (std::size_t index = 0; index < pairs.size(); ++index)
    {
        const auto &pair = pairs[index];
        values[index] = model == "sw" && pair.aqueous
            ? 0.0 : matrix[pair.first][pair.second];
    }
    return values;
}

Eos makeEos(const std::string &model, const std::array<double, 6> &parameters)
{
    Matrix matrix = baseMatrix(model);
    for (std::size_t index = 0; index < pairs.size(); ++index)
    {
        const auto &pair = pairs[index];
        if (!(model == "sw" && pair.aqueous))
        {
            matrix[pair.first][pair.second] = parameters[index];
            matrix[pair.second][pair.first] = parameters[index];
        }
    }

    MPMC::CompositionalMixture<Indices> mixture(
        CaseConfig::CommonFluid::criticalTemperature,
        CaseConfig::CommonFluid::criticalPressure,
        CaseConfig::CommonFluid::criticalVolume,
        CaseConfig::CommonFluid::acentricFactor,
        CaseConfig::CommonFluid::molarMass,
        matrix);
    Eos eos(
        CaseConfig::CommonFluid::eosOmegaA,
        CaseConfig::CommonFluid::eosOmegaB,
        std::move(mixture),
        model == "pr" ? CaseConfig::PrFluid::eosModelFlag :
            model == "sw" ? CaseConfig::SwFluid::eosModelFlag :
            CaseConfig::CpaFluid::eosModelFlag,
        CaseConfig::CommonFluid::eosU,
        CaseConfig::CommonFluid::eosW);

    if (model == "sw")
    {
        Eos::SoreideWhitsonOptions options;
        options.waterComponent = static_cast<int>(h2o);
        options.salinityMolality = 1.0;
        for (std::size_t component = 0; component < 5; ++component)
        {
            options.aqueousWaterBip[component] =
                [component, parameters](double temperature, double salinity) {
                    if (component == h2o)
                        return 0.0;
                    if (component == co2)
                    {
                        return MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(
                            temperature, salinity) + parameters[0];
                    }
                    const double baseline =
                        MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                            temperature,
                            CaseConfig::CommonFluid::criticalTemperature[component],
                            CaseConfig::CommonFluid::acentricFactor[component],
                            salinity);
                    return baseline + (component == ch4 ? parameters[1] : 0.0);
                };
        }
        eos.configureSoreideWhitson(std::move(options));
    }
    else if (model == "cpa")
    {
        Eos::CubicPlusAssociationOptions options;
        options.a0 = CaseConfig::CpaFluid::cpaA0;
        options.b = CaseConfig::CpaFluid::cpaB;
        options.c1 = CaseConfig::CpaFluid::cpaC1;
        options.associationEnergy = CaseConfig::CpaFluid::cpaAssociationEnergy;
        options.associationVolume = CaseConfig::CpaFluid::cpaAssociationVolume;
        options.donorSites = CaseConfig::CpaFluid::cpaDonorSites;
        options.acceptorSites = CaseConfig::CpaFluid::cpaAcceptorSites;
        options.crossAssociationEnergy = CaseConfig::CpaFluid::cpaCrossAssociationEnergy;
        options.crossAssociationVolume = CaseConfig::CpaFluid::cpaCrossAssociationVolume;
        options.radialDistribution = CaseConfig::CpaFluid::cpaRadialDistribution;
        eos.configureCubicPlusAssociation(std::move(options));
    }
    return eos;
}

std::vector<PointPrediction> predict(const Eos &eos, const Datum &datum)
{
    if (datum.kind == "density")
    {
        const auto gas = eos.phaseResult(
            datum.pressure, datum.temperature, datum.vapor,
            MPMC::CompositionalPhase::Gas, false);
        const double molarDensity = datum.pressure /
            (gas.compressibility * MPMC::units::gasConstant * datum.temperature);
        double molarMass = 0.0;
        for (std::size_t component = 0; component < 5; ++component)
            molarMass += datum.vapor[component] * CaseConfig::CommonFluid::molarMass[component];
        const double density = molarDensity * molarMass;
        return {{"mass_density_kg_m3", "mixture", datum.density, density,
            density / datum.density - 1.0}};
    }

    const auto liquidPhase = datum.pair.rfind("h2o_", 0) == 0
        ? MPMC::CompositionalPhase::Water : MPMC::CompositionalPhase::Oil;
    const auto liquid = eos.phaseResult(
        datum.pressure, datum.temperature, datum.liquid, liquidPhase, false);
    const auto vapor = eos.phaseResult(
        datum.pressure, datum.temperature, datum.vapor,
        MPMC::CompositionalPhase::Gas, false);

    const PairDefinition *definition = nullptr;
    for (const auto &pair : pairs)
        if (pair.name == datum.pair)
            definition = &pair;
    if (!definition)
        throw std::runtime_error("Unknown pair in common data: " + datum.pair);

    std::vector<std::size_t> active{definition->first, definition->second};
    if (datum.kind == "solubility")
        active = {co2};
    const std::array<std::string, 5> names{"H2O", "CO2", "CH4", "C2H6", "nC4H10"};
    std::vector<PointPrediction> result;
    for (const std::size_t component : active)
    {
        if (!(datum.liquid[component] > 0.0) || !(datum.vapor[component] > 0.0))
            continue;
        const double residual = std::log(
            liquid.fugacity[component] / vapor.fugacity[component]);
        result.push_back({"ln_fugacity_ratio", names[component], 0.0,
            residual, std::exp(residual) - 1.0});
    }
    return result;
}

double objectiveForPair(
    const std::string &model,
    const std::array<double, 6> &parameters,
    const std::vector<Datum> &data,
    const std::string &pair,
    const std::string &split)
{
    try
    {
        const Eos eos = makeEos(model, parameters);
        double sum = 0.0;
        std::size_t count = 0;
        for (const auto &datum : data)
        {
            if (datum.pair != pair || datum.split != split)
                continue;
            for (const auto &prediction : predict(eos, datum))
            {
                const double residual = prediction.quantity == "mass_density_kg_m3"
                    ? std::log(prediction.predicted / prediction.observed)
                    : prediction.predicted;
                sum += residual * residual;
                ++count;
            }
        }
        if (count == 0)
            throw std::runtime_error("No data for pair/split: " + pair + "/" + split);
        return sum / static_cast<double>(count);
    }
    catch (const std::exception &)
    {
        return 1.0e12;
    }
}

void writeParameters(
    const std::filesystem::path &path,
    const std::vector<FitRecord> &records)
{
    std::ofstream output(path);
    output << std::setprecision(17)
        << "model,pair,parameter_kind,initial,lower,upper,fitted,"
           "training_objective_initial,training_objective_final,evaluations\n";
    for (const auto &record : records)
    {
        output << record.model << ',' << record.pair << ',' << record.parameterKind << ','
               << record.initial << ',' << record.lower << ',' << record.upper << ','
               << record.fitted << ',' << record.trainingObjectiveInitial << ','
               << record.trainingObjectiveFinal << ',' << record.evaluations << '\n';
    }
}

void writePredictions(
    const std::filesystem::path &path,
    const std::vector<Datum> &data,
    const std::map<std::string, std::array<double, 6>> &fitted)
{
    std::ofstream output(path);
    output << std::setprecision(17)
        << "model,id,source,pair,kind,split,temperature_K,pressure_Pa,"
           "quantity,component,observed,predicted,relative_error\n";
    for (const auto &[model, parameters] : fitted)
    {
        const Eos eos = makeEos(model, parameters);
        for (const auto &datum : data)
        {
            try
            {
                for (const auto &prediction : predict(eos, datum))
                {
                    output << model << ',' << datum.id << ',' << datum.source << ','
                           << datum.pair << ',' << datum.kind << ',' << datum.split << ','
                           << datum.temperature << ',' << datum.pressure << ','
                           << prediction.quantity << ',' << prediction.component << ','
                           << prediction.observed << ',' << prediction.predicted << ','
                           << prediction.relativeError << '\n';
                }
            }
            catch (const std::exception &error)
            {
                std::cerr << "[common-fit] prediction failed: " << model << ' '
                          << datum.id << ": " << error.what() << '\n';
            }
        }
    }
}

void writeSummary(
    const std::filesystem::path &path,
    const std::vector<Datum> &data,
    const std::map<std::string, std::array<double, 6>> &fitted)
{
    std::ofstream output(path);
    output << std::setprecision(17)
        << "model,pair,split,points,predictions,rmse_log_ratio,"
           "mean_abs_percent_error,max_abs_percent_error\n";
    for (const auto &[model, parameters] : fitted)
    {
        const Eos eos = makeEos(model, parameters);
        for (const auto &pair : pairs)
        {
            for (const std::string split : {"training", "validation"})
            {
                std::vector<double> residuals;
                std::size_t pointCount = 0;
                for (const auto &datum : data)
                {
                    if (datum.pair != pair.name || datum.split != split)
                        continue;
                    ++pointCount;
                    try
                    {
                        for (const auto &prediction : predict(eos, datum))
                        {
                            residuals.push_back(
                                prediction.quantity == "mass_density_kg_m3"
                                ? std::log(prediction.predicted / prediction.observed)
                                : prediction.predicted);
                        }
                    }
                    catch (const std::exception &)
                    {
                        residuals.push_back(std::log(1.0e6));
                    }
                }
                double sumSquares = 0.0;
                double sumPercent = 0.0;
                double maxPercent = 0.0;
                for (const double residual : residuals)
                {
                    sumSquares += residual * residual;
                    const double percent = 100.0 * std::abs(std::exp(residual) - 1.0);
                    sumPercent += percent;
                    maxPercent = std::max(maxPercent, percent);
                }
                const double count = static_cast<double>(residuals.size());
                output << model << ',' << pair.name << ',' << split << ','
                       << pointCount << ',' << residuals.size() << ','
                       << std::sqrt(sumSquares / count) << ','
                       << sumPercent / count << ',' << maxPercent << '\n';
            }
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            std::cerr << "Usage: five_component_common_fit data.csv output_dir\n";
            return 2;
        }
        const auto data = readData(argv[1]);
        const std::filesystem::path outputDirectory(argv[2]);
        std::filesystem::create_directories(outputDirectory);

        std::vector<FitRecord> records;
        std::map<std::string, std::array<double, 6>> fitted;
        for (const std::string model : {"pr", "sw", "cpa"})
        {
            auto parameters = initialParameters(model);
            for (std::size_t pairIndex = 0; pairIndex < pairs.size(); ++pairIndex)
            {
                const auto &pair = pairs[pairIndex];
                const double initial = parameters[pairIndex];
                const MPMC::tools::RegressionParameter descriptor{
                    pair.name,
                    initial,
                    model == "sw" && pair.aqueous ? -0.20 : pair.lower,
                    model == "sw" && pair.aqueous ? 0.20 : pair.upper,
                    pair.step};
                MPMC::tools::PatternSearchOptions options;
                options.maximumIterations = 50;
                options.minimumRelativeStep = 1.0e-7;
                const auto result = MPMC::tools::boundedPatternSearch(
                    std::vector<MPMC::tools::RegressionParameter>{descriptor},
                    [&](const std::vector<double> &candidate) {
                        auto trial = parameters;
                        trial[pairIndex] = candidate[0];
                        return objectiveForPair(
                            model, trial, data, pair.name, "training");
                    },
                    options);
                parameters[pairIndex] = result.values[0];
                records.push_back({
                    model,
                    pair.name,
                    model == "sw" && pair.aqueous ? "aqueous_correlation_offset" : "kij",
                    initial,
                    descriptor.lowerBound,
                    descriptor.upperBound,
                    result.values[0],
                    result.initialObjective,
                    result.objective,
                    result.objectiveEvaluations});
                std::cout << "[common-fit] " << model << ' ' << pair.name
                          << " = " << std::setprecision(10) << result.values[0]
                          << ", objective " << result.initialObjective << " -> "
                          << result.objective << '\n';
            }
            fitted.emplace(model, parameters);
        }

        writeParameters(outputDirectory / "fitted_parameters.csv", records);
        writePredictions(outputDirectory / "point_predictions.csv", data, fitted);
        writeSummary(outputDirectory / "fit_summary.csv", data, fitted);
        std::cout << "[common-fit] wrote results to " << outputDirectory << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[common-fit] ERROR: " << error.what() << '\n';
        return 1;
    }
}
