/**
 * @file main.cpp
 * @brief H2O/角鲨烷零流动 PVT、BIP 与黏度标定。
 *
 * The executable reuses the production PR EOS, three-phase flash, volume
 * translation and LBC property model.  IAPWS-2008 is evaluated only as an
 * independent pure-water reference; it does not replace production physics.
 */
#include <common/math.hpp>
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/properties/aqueous_viscosity.hpp>
#include <natural/properties/compositional_properties.hpp>
#include <natural/thermo/aqueous_volume.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <tools/eos_parameter_regression.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
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
    2, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Mixture = MPMC::CompositionalMixture<Indices>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using PropertyModel = MPMC::CompositionalPropertyModel<Indices>;
using Composition = std::array<double, 2>;
using Matrix = std::array<std::array<double, 2>, 2>;

constexpr double referenceTemperature = 653.2;
constexpr double currentKij = 0.2395;
constexpr double currentInverseTemperatureSlope = -610.0;
constexpr double currentHeavyCriticalVolume = 3.05355324646748e-3;
constexpr double waterMolarMass = 0.01801528;
constexpr double squalaneMolarMass = 0.4228;
constexpr double missingValue = std::numeric_limits<double>::quiet_NaN();

constexpr std::array<double, 2> criticalTemperature{647.096, 855.9};
constexpr std::array<double, 2> criticalPressure{22.064e6, 0.7164e6};
constexpr std::array<double, 2> acentricFactor{0.3443, 1.255};
constexpr std::array<double, 2> molarMass{waterMolarMass, squalaneMolarMass};

struct LleDatum
{
    std::string id;
    std::string split;
    double temperature{};
    double pressure{};
    double xWaterRich{};
    double xSqualaneRich{};
    double accuracyWaterRich{};
    double accuracySqualaneRich{};
};

struct PvtDatum
{
    std::string id;
    std::string split;
    double temperature{};
    double pressure{};
    double density{};
    double viscosity{};
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

std::vector<LleDatum> readLle(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open LLE data: " + path.string());
    const auto column = readHeader(input);
    std::vector<LleDatum> rows;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto value = splitCsv(line);
        auto text = [&](const std::string &name) { return value.at(column.at(name)); };
        auto number = [&](const std::string &name) { return std::stod(text(name)); };
        rows.push_back({
            text("id"), text("role"), number("temperature_K"),
            number("pressure_MPa") * 1.0e6,
            number("x_water_rich"), number("x_water_squalane_rich"),
            number("reported_accuracy_water_rich"),
            number("reported_accuracy_squalane_rich")});
    }
    if (rows.size() != 4)
        throw std::runtime_error("Expected four squalane coexistence rows");
    return rows;
}

std::vector<PvtDatum> readPvt(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open PVT data: " + path.string());
    const auto column = readHeader(input);
    std::vector<PvtDatum> rows;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto value = splitCsv(line);
        auto text = [&](const std::string &name) { return value.at(column.at(name)); };
        auto number = [&](const std::string &name) { return std::stod(text(name)); };
        rows.push_back({
            text("id"), text("split"), number("temperature_K"),
            number("pressure_MPa") * 1.0e6, number("density_kg_m3"),
            number("viscosity_mPa_s") * 1.0e-3});
    }
    if (rows.size() != 24)
        throw std::runtime_error("Expected 24 squalane PVT/viscosity rows");
    return rows;
}

Mixture makeMixture(double heavyCriticalVolume, double constantKij = currentKij)
{
    const std::array<double, 2> criticalVolume{
        7.49587808839913e-5, heavyCriticalVolume};
    const Matrix bip{{{{0.0, constantKij}}, {{constantKij, 0.0}}}};
    return Mixture(
        criticalTemperature, criticalPressure, criticalVolume,
        acentricFactor, molarMass, bip);
}

Eos makeEos(
    double kReference,
    double inverseTemperatureSlope,
    double heavyVolumeTranslation,
    double heavyCriticalVolume)
{
    Eos eos(
        0.45724, 0.07780,
        makeMixture(heavyCriticalVolume, kReference),
        1, 1.0 + std::sqrt(2.0), 1.0 - std::sqrt(2.0), 1.0e-30);
    eos.configureBinaryInteractionFunction(
        [kReference, inverseTemperatureSlope](int i, int j, double temperature) {
            if (i == j)
                return 0.0;
            return kReference + inverseTemperatureSlope *
                (1.0 / temperature - 1.0 / referenceTemperature);
        });
    eos.configureVolumeTranslation({0.0, heavyVolumeTranslation});
    return eos;
}

double pureSqualaneDensity(
    double temperature,
    double pressure,
    double volumeTranslation,
    double heavyCriticalVolume)
{
    const Composition pure{0.0, 1.0};
    const auto eos = makeEos(
        currentKij, currentInverseTemperatureSlope,
        volumeTranslation, heavyCriticalVolume);
    const auto phase = eos.phaseResult(
        pressure, temperature, pure, MPMC::CompositionalPhase::Oil, false);
    return eos.molarDensity(
        pressure, temperature, pure, phase.compressibility,
        MPMC::CompositionalPhase::Oil) * squalaneMolarMass;
}

double pureSqualaneViscosity(
    double temperature,
    double pressure,
    double volumeTranslation,
    double heavyCriticalVolume)
{
    const Composition pure{0.0, 1.0};
    const auto mixture = makeMixture(heavyCriticalVolume);
    const PropertyModel properties(mixture);
    const auto eos = makeEos(
        currentKij, currentInverseTemperatureSlope,
        volumeTranslation, heavyCriticalVolume);
    const auto phase = eos.phaseResult(
        pressure, temperature, pure, MPMC::CompositionalPhase::Oil, false);
    const double molarDensity = eos.molarDensity(
        pressure, temperature, pure, phase.compressibility,
        MPMC::CompositionalPhase::Oil);
    return properties.viscosityFromMolarDensity(
        pure, molarDensity, temperature);
}

double pvtObjective(const std::vector<PvtDatum> &rows, double shift)
{
    double sum = 0.0;
    int count = 0;
    try
    {
        for (const auto &row : rows)
        {
            if (row.split != "training")
                continue;
            const double predicted = pureSqualaneDensity(
                row.temperature, row.pressure, shift,
                currentHeavyCriticalVolume);
            const double residual = std::log(predicted / row.density);
            sum += residual * residual;
            ++count;
        }
    }
    catch (const std::exception &)
    {
        return 1.0e12;
    }
    return count > 0 ? sum / count : 1.0e12;
}

double viscosityObjective(
    const std::vector<PvtDatum> &rows,
    double volumeTranslation,
    double log10CriticalVolume)
{
    const double criticalVolume = std::pow(10.0, log10CriticalVolume);
    double sum = 0.0;
    int count = 0;
    try
    {
        for (const auto &row : rows)
        {
            if (row.split != "training")
                continue;
            const double predicted = pureSqualaneViscosity(
                row.temperature, row.pressure,
                volumeTranslation, criticalVolume);
            const double residual = std::log(predicted / row.viscosity);
            sum += residual * residual;
            ++count;
        }
    }
    catch (const std::exception &)
    {
        return 1.0e12;
    }
    return count > 0 ? sum / count : 1.0e12;
}

std::array<double, 2> endpointFugacityResidual(
    const Eos &eos,
    const LleDatum &row)
{
    const Composition waterRich{row.xWaterRich, 1.0 - row.xWaterRich};
    const Composition squalaneRich{
        row.xSqualaneRich, 1.0 - row.xSqualaneRich};
    const auto first = eos.phaseResult(
        row.pressure, row.temperature, waterRich, false, true);
    const auto second = eos.phaseResult(
        row.pressure, row.temperature, squalaneRich, true, true);
    return {
        std::log(first.fugacity[0] / second.fugacity[0]),
        std::log(first.fugacity[1] / second.fugacity[1])};
}

double singleStateFlashObjective(
    const LleDatum &row,
    double kij)
{
    try
    {
        const auto eos = makeEos(
            kij, 0.0, 0.0, currentHeavyCriticalVolume);
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = 0;
        const Flash flash(eos, options);
        const auto result = flash.flash(
            row.pressure, row.temperature, Composition{0.96, 0.04});
        std::vector<double> waterFractions;
        if (result.converged)
        {
            for (std::size_t phase = 0; phase < 3; ++phase)
            {
                if (result.phaseMoleFraction[phase] > options.phaseFractionTolerance)
                    waterFractions.push_back(result.composition[phase][0]);
            }
        }
        if (waterFractions.size() != 2)
            return 1.0e12;
        std::sort(
            waterFractions.begin(), waterFractions.end(),
            std::greater<double>());
        const double waterResidual =
            (waterFractions[0] - row.xWaterRich) / row.accuracyWaterRich;
        const double squalaneResidual =
            (waterFractions[1] - row.xSqualaneRich)
            / row.accuracySqualaneRich;
        return 0.5 * (waterResidual * waterResidual
                      + squalaneResidual * squalaneResidual);
    }
    catch (const std::exception &)
    {
        return 1.0e12;
    }
}

double bipObjective(
    const std::vector<LleDatum> &rows,
    double kReference,
    double inverseTemperatureSlope)
{
    double sum = 0.0;
    int count = 0;
    for (const auto &row : rows)
    {
        if (row.split != "calibration")
            continue;
        const double kij = kReference + inverseTemperatureSlope *
            (1.0 / row.temperature - 1.0 / referenceTemperature);
        if (kij < -0.2 || kij > 0.8)
            return 1.0e12;
        const double value = singleStateFlashObjective(row, kij);
        if (!(value < 1.0e12))
            return 1.0e12;
        sum += value;
        ++count;
    }
    return count > 0 ? sum / count : 1.0e12;
}

double iapws2008ViscosityIndustrial(double temperature, double density)
{
    const double reducedTemperature = temperature / 647.096;
    const double reducedDensity = density / 322.0;
    constexpr std::array<double, 4> h{
        1.67752, 2.20462, 0.6366564, -0.241605};
    double denominator = 0.0;
    for (std::size_t i = 0; i < h.size(); ++i)
        denominator += h[i] / std::pow(reducedTemperature, static_cast<int>(i));
    const double mu0 = 100.0 * std::sqrt(reducedTemperature) / denominator;

    struct Coefficient { int i; int j; double value; };
    constexpr std::array<Coefficient, 21> coefficients{{
        {0, 0, 5.20094e-1}, {1, 0, 8.50895e-2}, {2, 0, -1.08374},
        {3, 0, -2.89555e-1}, {0, 1, 2.22531e-1}, {1, 1, 9.99115e-1},
        {2, 1, 1.88797}, {3, 1, 1.26613}, {5, 1, 1.20573e-1},
        {0, 2, -2.81378e-1}, {1, 2, -9.06851e-1}, {2, 2, -7.72479e-1},
        {3, 2, -4.89837e-1}, {4, 2, -2.57040e-1}, {0, 3, 1.61913e-1},
        {1, 3, 2.57399e-1}, {0, 4, -3.25372e-2}, {3, 4, 6.98452e-2},
        {4, 5, 8.72102e-3}, {3, 6, -4.35673e-3}, {5, 6, -5.93264e-4}}};
    const double temperatureVariable = 1.0 / reducedTemperature - 1.0;
    const double densityVariable = reducedDensity - 1.0;
    double exponentSeries = 0.0;
    for (const auto &coefficient : coefficients)
    {
        exponentSeries += coefficient.value
            * std::pow(temperatureVariable, coefficient.i)
            * std::pow(densityVariable, coefficient.j);
    }
    const double mu1 = std::exp(reducedDensity * exponentSeries);
    return mu0 * mu1 * 1.0e-6;
}

double averageAbsoluteRelativePercent(
    const std::vector<PvtDatum> &rows,
    const std::string &split,
    const std::function<double(const PvtDatum &)> &predict,
    bool viscosity)
{
    double sum = 0.0;
    int count = 0;
    for (const auto &row : rows)
    {
        if (row.split != split)
            continue;
        const double target = viscosity ? row.viscosity : row.density;
        sum += std::abs(predict(row) / target - 1.0);
        ++count;
    }
    return count > 0 ? 100.0 * sum / count : missingValue;
}

void writeFitParameters(
    const std::filesystem::path &path,
    const MPMC::tools::RegressionResult &pvtFit,
    const MPMC::tools::RegressionResult &viscosityFit,
    const MPMC::tools::RegressionResult &bipFit,
    double initialKReference,
    double initialInverseSlope)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create fit parameter output");
    out << std::setprecision(17)
        << "parameter,initial,lower,upper,fitted,objective_initial,objective_final,evaluations\n"
        << "squalane_volume_translation_m3_mol,0,-0.0005,0.0005,"
        << pvtFit.values.at(0) << ',' << pvtFit.initialObjective << ','
        << pvtFit.objective << ',' << pvtFit.objectiveEvaluations << '\n'
        << "log10_squalane_effective_critical_volume,-2.515194,-5,-2,"
        << viscosityFit.values.at(0) << ',' << viscosityFit.initialObjective << ','
        << viscosityFit.objective << ',' << viscosityFit.objectiveEvaluations << '\n'
        << "squalane_effective_critical_volume_m3_mol,0.00305355324646748,0.00001,0.01,"
        << std::pow(10.0, viscosityFit.values.at(0)) << ','
        << viscosityFit.initialObjective << ',' << viscosityFit.objective << ','
        << viscosityFit.objectiveEvaluations << '\n'
        << "k_H2O_squalane_at_653.2K," << initialKReference << ",-0.2,0.8,"
        << bipFit.values.at(0) << ',' << bipFit.initialObjective << ','
        << bipFit.objective << ',' << bipFit.objectiveEvaluations << '\n'
        << "b_inverse_temperature_K," << initialInverseSlope << ",-5000,5000,"
        << bipFit.values.at(1) << ',' << bipFit.initialObjective << ','
        << bipFit.objective << ',' << bipFit.objectiveEvaluations << '\n';
}

void writePvtPredictions(
    const std::filesystem::path &path,
    const std::vector<PvtDatum> &rows,
    double shift,
    double fittedCriticalVolume)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create PVT prediction output");
    out << std::setprecision(17)
        << "id,split,temperature_K,pressure_MPa,target_density_kg_m3,"
           "unshifted_PR_density_kg_m3,translated_PR_density_kg_m3,"
           "target_viscosity_Pa_s,current_LBC_viscosity_Pa_s,"
           "calibrated_LBC_viscosity_Pa_s\n";
    for (const auto &row : rows)
    {
        out << row.id << ',' << row.split << ',' << row.temperature << ','
            << row.pressure / 1.0e6 << ',' << row.density << ','
            << pureSqualaneDensity(
                row.temperature, row.pressure, 0.0,
                currentHeavyCriticalVolume) << ','
            << pureSqualaneDensity(
                row.temperature, row.pressure, shift,
                currentHeavyCriticalVolume) << ','
            << row.viscosity << ','
            << pureSqualaneViscosity(
                row.temperature, row.pressure, 0.0,
                currentHeavyCriticalVolume) << ','
            << pureSqualaneViscosity(
                row.temperature, row.pressure, shift,
                fittedCriticalVolume) << '\n';
    }
}

void writeLlePredictions(
    const std::filesystem::path &path,
    const std::vector<LleDatum> &rows,
    double kReference,
    double inverseTemperatureSlope,
    double shift,
    double fittedCriticalVolume)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create LLE prediction output");
    out << std::setprecision(17)
        << "variant,id,split,temperature_K,pressure_MPa,kij,target_x_water_rich,"
           "target_x_squalane_rich,predicted_x_water_rich,"
           "predicted_x_squalane_rich,phase_count,converged,"
           "endpoint_logf_residual_h2o,endpoint_logf_residual_squalane,"
           "water_rich_density_kg_m3,squalane_rich_density_kg_m3,"
           "water_rich_viscosity_Pa_s,squalane_rich_viscosity_Pa_s\n";
    const auto emitVariant = [&](
        const std::string &variant,
        double variantKReference,
        double variantSlope,
        double variantShift,
        double variantCriticalVolume) {
        const auto eos = makeEos(
            variantKReference, variantSlope,
            variantShift, variantCriticalVolume);
        const PropertyModel properties(makeMixture(variantCriticalVolume));
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = 0;
        const Flash flash(eos, options);
        for (const auto &row : rows)
        {
            const auto endpointResidual = endpointFugacityResidual(eos, row);
            const auto result = flash.flash(
                row.pressure, row.temperature, Composition{0.96, 0.04});
            std::vector<std::size_t> active;
            for (std::size_t phase = 0; phase < 3; ++phase)
            {
                if (result.phaseMoleFraction[phase] > options.phaseFractionTolerance)
                    active.push_back(phase);
            }
            std::sort(active.begin(), active.end(), [&](std::size_t a, std::size_t b) {
                return result.composition[a][0] > result.composition[b][0];
            });
            double xWaterRich = missingValue;
            double xSqualaneRich = missingValue;
            double densityWaterRich = missingValue;
            double densitySqualaneRich = missingValue;
            double viscosityWaterRich = missingValue;
            double viscositySqualaneRich = missingValue;
            if (active.size() == 2)
            {
                const auto evaluate = [&](std::size_t slot) {
                    const auto role = static_cast<MPMC::CompositionalPhase>(slot);
                    const double molarDensity = eos.molarDensity(
                        row.pressure, row.temperature, result.composition[slot],
                        result.compressibility[slot], role);
                    const double mixtureMass =
                        result.composition[slot][0] * waterMolarMass
                        + result.composition[slot][1] * squalaneMolarMass;
                    return std::array<double, 2>{
                        molarDensity * mixtureMass,
                        properties.viscosityFromMolarDensity(
                            result.composition[slot], molarDensity, row.temperature)};
                };
                xWaterRich = result.composition[active[0]][0];
                xSqualaneRich = result.composition[active[1]][0];
                const auto waterValues = evaluate(active[0]);
                const auto squalaneValues = evaluate(active[1]);
                densityWaterRich = waterValues[0];
                viscosityWaterRich = waterValues[1];
                densitySqualaneRich = squalaneValues[0];
                viscositySqualaneRich = squalaneValues[1];
            }
            out << variant << ',' << row.id << ',' << row.split << ','
                << row.temperature << ',' << row.pressure / 1.0e6 << ','
                << eos.binaryInteractionCoefficient(0, 1, row.temperature) << ','
                << row.xWaterRich << ',' << row.xSqualaneRich << ','
                << xWaterRich << ',' << xSqualaneRich << ',' << active.size() << ','
                << (result.converged ? 1 : 0) << ',' << endpointResidual[0] << ','
                << endpointResidual[1] << ',' << densityWaterRich << ','
                << densitySqualaneRich << ',' << viscosityWaterRich << ','
                << viscositySqualaneRich << '\n';
        }
    };
    emitVariant(
        "current_case", currentKij, currentInverseTemperatureSlope,
        0.0, currentHeavyCriticalVolume);
    emitVariant(
        "flash_composition_fit", kReference, inverseTemperatureSlope,
        shift, fittedCriticalVolume);
}

void writeWaterCheck(
    const std::filesystem::path &path,
    const std::vector<LleDatum> &rows)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create water viscosity output");
    out << std::setprecision(17)
        << "id,temperature_K,pressure_MPa,if97_density_kg_m3,"
           "iapws2008_viscosity_Pa_s,mcbride_wright_viscosity_Pa_s,"
           "mcbride_relative_error_percent\n";
    const MPMC::McBrideWright2015AqueousViscosity<2> mcBride(0, 1, 0.0);
    const Composition pureWater{1.0, 0.0};
    for (const auto &row : rows)
    {
        const double density = MPMC::IapwsIf97WaterDensity::density(
            row.pressure, row.temperature);
        const double reference = iapws2008ViscosityIndustrial(
            row.temperature, density);
        const double current = mcBride.viscosity(
            row.pressure, row.temperature, pureWater);
        out << row.id << ',' << row.temperature << ','
            << row.pressure / 1.0e6 << ',' << density << ',' << reference << ','
            << current << ',' << 100.0 * (current / reference - 1.0) << '\n';
    }
}

void writeSummary(
    const std::filesystem::path &path,
    const std::vector<PvtDatum> &pvtRows,
    const std::vector<LleDatum> &lleRows,
    double shift,
    double fittedCriticalVolume,
    double kReference,
    double inverseTemperatureSlope)
{
    const auto densityPredict = [&](const PvtDatum &row) {
        return pureSqualaneDensity(
            row.temperature, row.pressure, shift,
            currentHeavyCriticalVolume);
    };
    const auto viscosityPredict = [&](const PvtDatum &row) {
        return pureSqualaneViscosity(
            row.temperature, row.pressure, shift,
            fittedCriticalVolume);
    };
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create summary output");
    out << std::setprecision(17) << "metric,value\n"
        << "squalane_volume_translation_m3_mol," << shift << '\n'
        << "squalane_effective_critical_volume_m3_mol,"
        << fittedCriticalVolume << '\n'
        << "density_training_AARD_percent,"
        << averageAbsoluteRelativePercent(
            pvtRows, "training", densityPredict, false) << '\n'
        << "density_validation_AARD_percent,"
        << averageAbsoluteRelativePercent(
            pvtRows, "validation", densityPredict, false) << '\n'
        << "viscosity_training_AARD_percent,"
        << averageAbsoluteRelativePercent(
            pvtRows, "training", viscosityPredict, true) << '\n'
        << "viscosity_validation_AARD_percent,"
        << averageAbsoluteRelativePercent(
            pvtRows, "validation", viscosityPredict, true) << '\n'
        << "k_ref_653.2K," << kReference << '\n'
        << "b_inverse_temperature_K," << inverseTemperatureSlope << '\n';

    const auto eos = makeEos(
        kReference, inverseTemperatureSlope, 0.0,
        currentHeavyCriticalVolume);
    for (const auto &row : lleRows)
    {
        const auto residual = endpointFugacityResidual(eos, row);
        out << "endpoint_logf_rms_" << row.id << ','
            << std::sqrt(0.5 * (residual[0] * residual[0]
                                + residual[1] * residual[1])) << '\n';
    }
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 4)
        {
            std::cerr << "Usage: scw_binary_calibration LLE.csv PVT.csv OUTPUT_DIR\n";
            return 2;
        }
        const auto lleRows = readLle(argv[1]);
        const auto pvtRows = readPvt(argv[2]);
        const double waterVerification = iapws2008ViscosityIndustrial(298.15, 998.0);
        if (std::abs(waterVerification / 889.735100e-6 - 1.0) > 1.0e-9)
            throw std::runtime_error("IAPWS-2008 viscosity verification point failed");
        const std::filesystem::path output(argv[3]);
        if (std::filesystem::exists(output))
            throw std::runtime_error("Output directory already exists: " + output.string());
        std::filesystem::create_directories(output);

        MPMC::tools::PatternSearchOptions options;
        options.maximumIterations = 180;
        options.minimumRelativeStep = 1.0e-10;
        options.improvementTolerance = 1.0e-15;

        const auto pvtFit = MPMC::tools::boundedPatternSearch(
            {{"squalane_volume_translation_m3_mol", 0.0, -5.0e-4, 5.0e-4, 1.0e-5}},
            [&](const std::vector<double> &x) {
                return pvtObjective(pvtRows, x[0]);
            }, options);
        const double volumeTranslation = pvtFit.values.at(0);

        const auto viscosityFit = MPMC::tools::boundedPatternSearch(
            {{"log10_squalane_effective_critical_volume",
              std::log10(currentHeavyCriticalVolume), -5.0, -2.0, 0.1}},
            [&](const std::vector<double> &x) {
                return viscosityObjective(pvtRows, volumeTranslation, x[0]);
            }, options);
        const double fittedCriticalVolume = std::pow(10.0, viscosityFit.values.at(0));

        std::array<double, 2> independentlyFittedK{
            missingValue, missingValue};
        std::size_t calibrationIndex = 0;
        for (const auto &row : lleRows)
        {
            if (row.split != "calibration")
                continue;
            double initialK = currentKij;
            double best = std::numeric_limits<double>::infinity();
            for (int index = 0; index <= 100; ++index)
            {
                const double candidate = -0.2 + 0.01 * index;
                const double value = singleStateFlashObjective(row, candidate);
                if (value < best)
                {
                    best = value;
                    initialK = candidate;
                }
            }
            const auto stateFit = MPMC::tools::boundedPatternSearch(
                {{"k_H2O_squalane", initialK, -0.2, 0.8, 0.01}},
                [&](const std::vector<double> &x) {
                    return singleStateFlashObjective(row, x[0]);
                }, options);
            if (calibrationIndex >= independentlyFittedK.size())
                throw std::runtime_error("Unexpected number of calibration temperatures");
            independentlyFittedK[calibrationIndex++] = stateFit.values.at(0);
        }
        if (calibrationIndex != 2)
            throw std::runtime_error("Expected two calibration coexistence rows");
        const LleDatum *lowTemperatureCalibration = nullptr;
        const LleDatum *referenceCalibration = nullptr;
        for (const auto &row : lleRows)
        {
            if (row.split != "calibration")
                continue;
            if (std::abs(row.temperature - referenceTemperature) < 1.0e-9)
                referenceCalibration = &row;
            else
                lowTemperatureCalibration = &row;
        }
        if (!lowTemperatureCalibration || !referenceCalibration)
            throw std::runtime_error("Calibration temperatures are not the expected 637.2/653.2 K pair");
        // Rows are sorted by temperature in the source table.
        const double initialKReference = independentlyFittedK[1];
        const double initialInverseSlope =
            (independentlyFittedK[0] - initialKReference)
            / (1.0 / lowTemperatureCalibration->temperature
               - 1.0 / referenceCalibration->temperature);
        if (initialInverseSlope < -5000.0 || initialInverseSlope > 5000.0)
            throw std::runtime_error("Independent flash fits imply a BIP slope outside the declared bounds");

        const auto bipFit = MPMC::tools::boundedPatternSearch(
            {
                {"k_H2O_squalane_at_653.2K", initialKReference,
                 -0.2, 0.8, 0.002},
                {"b_inverse_temperature_K", initialInverseSlope,
                 -5000.0, 5000.0, 20.0}
            },
            [&](const std::vector<double> &x) {
                return bipObjective(lleRows, x[0], x[1]);
            }, options);

        writeFitParameters(
            output / "fit_parameters.csv", pvtFit, viscosityFit, bipFit,
            initialKReference, initialInverseSlope);
        writePvtPredictions(
            output / "squalane_pvt_viscosity_predictions.csv", pvtRows,
            volumeTranslation, fittedCriticalVolume);
        writeLlePredictions(
            output / "squalane_lle_predictions.csv", lleRows,
            bipFit.values.at(0), bipFit.values.at(1),
            volumeTranslation, fittedCriticalVolume);
        writeWaterCheck(
            output / "water_viscosity_check.csv", lleRows);
        writeSummary(
            output / "summary.csv", pvtRows, lleRows,
            volumeTranslation, fittedCriticalVolume,
            bipFit.values.at(0), bipFit.values.at(1));

        std::cout << std::setprecision(12)
            << "volume_translation=" << volumeTranslation << " m3/mol\n"
            << "effective_critical_volume=" << fittedCriticalVolume << " m3/mol\n"
            << "k_ref=" << bipFit.values.at(0)
            << " b_invT=" << bipFit.values.at(1) << " K\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
