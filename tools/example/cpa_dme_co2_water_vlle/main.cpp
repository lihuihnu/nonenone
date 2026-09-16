/**
 * @file main.cpp
 * @brief 基于 Laursen 等公开数据的 CO2/DME/H2O 三相 CPA 基准程序。
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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using Config = MPMC::CompositionalModelConfig<
    3, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, 3>;

constexpr std::size_t water = 0;
constexpr std::size_t dme = 1;
constexpr double temperature = 308.15;

struct ExperimentalPoint
{
    double pressureBar{0.0};
    Composition lower{}; // water-rich liquid, order H2O/DME/CO2
    Composition upper{}; // DME-rich liquid
    Composition vapour{};
};

std::vector<std::string> split(const std::string &line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
    {
        if (!field.empty() && field.back() == '\r')
            field.pop_back();
        fields.push_back(field);
    }
    return fields;
}

void normalize(Composition &x)
{
    const double sum = x[0] + x[1] + x[2];
    if (!(sum > 0.0))
        throw std::runtime_error("Cannot normalize an empty composition.");
    for (double &value : x)
        value /= sum;
}

std::vector<ExperimentalPoint> readReference(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open reference CSV: " + path.string());
    std::string line;
    std::getline(input, line);
    std::vector<ExperimentalPoint> points;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto fields = split(line);
        if (fields.size() != 11)
            throw std::runtime_error("Malformed reference row: " + line);
        const double rowTemperature = std::stod(fields[0]);
        if (std::abs(rowTemperature - temperature) > 1.0e-10)
            throw std::runtime_error("This benchmark expects the 308.15 K isotherm.");
        ExperimentalPoint point;
        point.pressureBar = std::stod(fields[1]);
        point.lower = {
            std::stod(fields[4]), std::stod(fields[3]), std::stod(fields[2])};
        point.upper = {
            std::stod(fields[7]), std::stod(fields[6]), std::stod(fields[5])};
        point.vapour = {
            std::stod(fields[10]), std::stod(fields[9]), std::stod(fields[8])};
        normalize(point.lower);
        normalize(point.upper);
        normalize(point.vapour);
        points.push_back(point);
    }
    return points;
}

Eos makeCpa()
{
    // Component order: H2O / DME / CO2. DME critical metadata only supplies
    // Wilson initialization; its fitted CPA a0/b/c1 below control properties.
    const std::array<double, 3> tc{647.29, 400.1, 304.1282};
    const std::array<double, 3> pc{22.064e6, 5.3368e6, 7.3773e6};
    const std::array<double, 3> vc{5.60e-5, 1.70e-4, 9.40e-5};
    const std::array<double, 3> omega{0.3443, 0.2000, 0.22394};
    const std::array<double, 3> mw{0.01801528, 0.04606844, 0.0440095};
    const std::array<std::array<double, 3>, 3> kij{{
        {{0.0, -0.160, -0.066}},
        {{-0.160, 0.0, -0.016}},
        {{-0.066, -0.016, 0.0}}
    }};
    Eos eos(
        0.42748, 0.08664,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        1, 1.0, 0.0, 1.0e-30);

    Eos::CubicPlusAssociationOptions cpa;
    constexpr double R = MPMC::units::gasConstant;
    for (std::size_t i = 0; i < 3; ++i)
    {
        cpa.a0[i] = 0.42748 * R * R * tc[i] * tc[i] / pc[i];
        cpa.b[i] = 0.08664 * R * tc[i] / pc[i];
        cpa.c1[i] = 0.480 + 1.574 * omega[i] - 0.176 * omega[i] * omega[i];
    }

    // Folas, Table 1.4: standard 4C water.
    cpa.a0[water] = 0.12277;
    cpa.b[water] = 1.4515e-5;
    cpa.c1[water] = 0.67359;
    cpa.associationEnergy[water] = 16655.0;
    cpa.associationVolume[water] = 0.0692;
    cpa.donorSites[water] = 2;
    cpa.acceptorSites[water] = 2;

    // Folas, Table 7.1: inert DME pure CPA parameters. Unit conversion:
    // 8.4354 bar L^2/mol^2 = 0.84354 Pa m^6/mol^2.
    cpa.a0[dme] = 0.84354;
    cpa.b[dme] = 4.96e-5;
    cpa.c1[dme] = 0.72125;

    // DME is not self-associating, but it has one acceptor site for water-DME
    // solvation. Folas mCR-1: epsilon_cross=epsilon_water/2 and
    // beta_cross=BETCR=0.2877. The same values are placed in DME's pure-site
    // storage solely to satisfy the generic parameter-integrity gate. They
    // cannot create DME-DME bonds because DME has zero donor sites, and the
    // explicit water->DME matrix entries below control the only active pair.
    cpa.acceptorSites[dme] = 1;
    cpa.associationEnergy[dme] = 0.5 * 16655.0;
    cpa.associationVolume[dme] = 0.2877;
    cpa.crossAssociationEnergy[water][dme] = 0.5 * 16655.0;
    cpa.crossAssociationVolume[water][dme] = 0.2877;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

Composition equalPhaseFeed(const ExperimentalPoint &point)
{
    Composition z{};
    for (std::size_t i = 0; i < 3; ++i)
        z[i] = (point.lower[i] + point.upper[i] + point.vapour[i]) / 3.0;
    normalize(z);
    return z;
}

double maxMassClosure(
    const Composition &z,
    const Flash::Result &result)
{
    double maximum = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
    {
        double reconstructed = 0.0;
        for (std::size_t phase = 0; phase < 3; ++phase)
        {
            reconstructed += result.phaseMoleFraction[phase] *
                result.composition[phase][component];
        }
        maximum = std::max(maximum, std::abs(reconstructed - z[component]));
    }
    return maximum;
}

double maxLogFugacitySpread(
    const Eos &eos,
    double pressure,
    const Flash::Result &result)
{
    const auto oil = eos.phaseResult(
        pressure, temperature, result.composition[0],
        MPMC::CompositionalPhase::Oil, false);
    const auto gas = eos.phaseResult(
        pressure, temperature, result.composition[1],
        MPMC::CompositionalPhase::Gas, false);
    const auto aqueous = eos.phaseResult(
        pressure, temperature, result.composition[2],
        MPMC::CompositionalPhase::Water, false);
    double maximum = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
    {
        const std::array<double, 3> logF{
            std::log(oil.fugacity[component]),
            std::log(gas.fugacity[component]),
            std::log(aqueous.fugacity[component])};
        maximum = std::max(maximum,
            *std::max_element(logF.begin(), logF.end()) -
            *std::min_element(logF.begin(), logF.end()));
    }
    return maximum;
}

void writePhaseRows(
    std::ofstream &output,
    double pressureBar,
    const std::string &phase,
    const Composition &experiment,
    const Composition &calculated)
{
    static const std::array<std::string, 3> names{"H2O", "DME", "CO2"};
    for (std::size_t component = 0; component < 3; ++component)
    {
        output << pressureBar << ',' << phase << ',' << names[component] << ','
               << experiment[component] << ',' << calculated[component] << ','
               << std::abs(calculated[component] - experiment[component]) << '\n';
    }
}

int run(const std::filesystem::path &reference, const std::filesystem::path &outputDir)
{
    const auto points = readReference(reference);
    std::filesystem::create_directories(outputDir);
    const Eos eos = makeCpa();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = static_cast<int>(water);
    options.maximumIterations = 180;
    options.maximumStabilityIterations = 120;
    const Flash flash(eos, options);

    std::ofstream summary(outputDir / "flash_results.csv");
    std::ofstream comparison(outputDir / "phase_comparison.csv");
    if (!summary || !comparison)
        throw std::runtime_error("Cannot create benchmark output CSV files.");
    summary << std::scientific << std::setprecision(12);
    comparison << std::scientific << std::setprecision(12);
    summary << "temperature_K,pressure_bar,converged,phase_code,iterations,"
               "z_H2O,z_DME,z_CO2,beta_oil,beta_gas,beta_water,"
               "max_mass_closure,max_log_fugacity_spread\n";
    comparison << "pressure_bar,phase,component,experimental,calculated,abs_error\n";

    std::size_t failures = 0;
    double maximumCompositionError = 0.0;
    double sumAbsoluteError = 0.0;
    std::size_t errorCount = 0;
    std::array<std::array<double, 3>, 3> relativeErrorSum{};
    std::array<std::array<std::size_t, 3>, 3> relativeErrorCount{};
    for (const auto &point : points)
    {
        const Composition z = equalPhaseFeed(point);
        const double pressure = point.pressureBar * 1.0e5;
        const auto result = flash.flash(pressure, temperature, z);
        const bool threePhase = result.converged &&
            result.presence.bits() == MPMC::PhasePresence::allBits;
        failures += threePhase ? 0u : 1u;
        const double closure = result.converged ? maxMassClosure(z, result) :
            std::numeric_limits<double>::quiet_NaN();
        const double fugacity = threePhase ? maxLogFugacitySpread(eos, pressure, result) :
            std::numeric_limits<double>::quiet_NaN();
        summary << temperature << ',' << point.pressureBar << ','
                << (result.converged ? 1 : 0) << ','
                << (result.converged ? result.presence.bits() : 0) << ','
                << result.iterations << ',' << z[0] << ',' << z[1] << ',' << z[2] << ','
                << result.phaseMoleFraction[0] << ',' << result.phaseMoleFraction[1] << ','
                << result.phaseMoleFraction[2] << ',' << closure << ',' << fugacity << '\n';

        if (!threePhase)
        {
            std::cout << point.pressureBar << " bar: FAIL (phase code "
                      << (result.converged ? result.presence.bits() : 0) << ")\n";
            continue;
        }

        // Public phase slots: oil=DME-rich upper liquid, gas=vapour,
        // water=water-rich lower liquid.
        writePhaseRows(comparison, point.pressureBar, "water-rich",
                       point.lower, result.composition[2]);
        writePhaseRows(comparison, point.pressureBar, "DME-rich",
                       point.upper, result.composition[0]);
        writePhaseRows(comparison, point.pressureBar, "vapor",
                       point.vapour, result.composition[1]);
        const std::array<Composition, 3> experiment{
            point.upper, point.vapour, point.lower};
        for (std::size_t phase = 0; phase < 3; ++phase)
        {
            for (std::size_t component = 0; component < 3; ++component)
            {
                const double error = std::abs(
                    result.composition[phase][component] - experiment[phase][component]);
                maximumCompositionError = std::max(maximumCompositionError, error);
                sumAbsoluteError += error;
                ++errorCount;
                if (experiment[phase][component] > 0.0)
                {
                    relativeErrorSum[phase][component] +=
                        100.0 * error / experiment[phase][component];
                    ++relativeErrorCount[phase][component];
                }
            }
        }
        std::cout << std::fixed << std::setprecision(1) << point.pressureBar
                  << " bar: three-phase PASS, beta="
                  << std::setprecision(6) << result.phaseMoleFraction[0] << '/'
                  << result.phaseMoleFraction[1] << '/'
                  << result.phaseMoleFraction[2] << ", closure="
                  << std::scientific << closure << ", fugacity=" << fugacity << '\n';
    }

    const double mae = errorCount > 0 ?
        sumAbsoluteError / static_cast<double>(errorCount) :
        std::numeric_limits<double>::quiet_NaN();
    std::ofstream metrics(outputDir / "metrics.csv");
    metrics << "points,three_phase_pass,failures,composition_mae,max_abs_composition_error\n"
            << points.size() << ',' << (points.size() - failures) << ',' << failures << ','
            << std::scientific << std::setprecision(12) << mae << ','
            << maximumCompositionError << '\n';

    std::ofstream aard(outputDir / "aard_by_phase_component.csv");
    if (!aard)
        throw std::runtime_error("Cannot create AARD output CSV file.");
    aard << "phase,component,points,our_aard_percent\n";
    const std::array<std::string, 3> phaseNames{"DME-rich", "vapor", "water-rich"};
    const std::array<std::string, 3> componentNames{"H2O", "DME", "CO2"};
    aard << std::scientific << std::setprecision(12);
    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        for (std::size_t component = 0; component < 3; ++component)
        {
            const auto count = relativeErrorCount[phase][component];
            const double value = count > 0 ? relativeErrorSum[phase][component] /
                static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
            aard << phaseNames[phase] << ',' << componentNames[component] << ','
                 << count << ',' << value << '\n';
        }
    }
    std::cout << "Composition MAE=" << mae
              << ", max abs error=" << maximumCompositionError << '\n';
    return failures == 0 ? 0 : 2;
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            std::cerr << "usage: cpa_dme_co2_water_vlle REFERENCE_CSV OUTPUT_DIR\n";
            return 1;
        }
        return run(argv[1], argv[2]);
    }
    catch (const std::exception &error)
    {
        std::cerr << "CPA DME VLLE benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
