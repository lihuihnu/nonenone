/**
 * @file main.cpp
 * @brief 公开数据与独立软件三 EOS/Flash 基准的生产内核批处理入口。
 */
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/soreide_whitson.hpp>
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
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{

using OgConfig = MPMC::CompositionalModelConfig<2, false, false>;
using OgIndices = MPMC::ScalarIndices<OgConfig>;
using OgEos = MPMC::CubicEquationOfState<OgIndices>;
using OgFlash = MPMC::CubicThreePhaseFlash<OgIndices>;
using OgComposition = std::array<double, 2>;

using SwConfig = MPMC::CompositionalModelConfig<
    2, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using SwIndices = MPMC::ScalarIndices<SwConfig>;
using SwEos = MPMC::CubicEquationOfState<SwIndices>;
using SwFlash = MPMC::CubicThreePhaseFlash<SwIndices>;
using SwComposition = std::array<double, 2>;

struct Request
{
    std::string id;
    std::string model;
    std::string operation;
    double temperature{0.0};
    double pressure{0.0};
    double salinity{0.0};
    double z0{0.5};
    double z1{0.5};
};

struct Output
{
    std::string id;
    bool converged{false};
    int phaseCode{0};
    int iterations{0};
    std::array<double, 3> beta{};
    std::array<std::array<double, 2>, 3> x{};
    std::array<double, 3> zFactor{};
    std::array<double, 3> molarDensity{};
    std::array<double, 2> lnPhi{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    double value{std::numeric_limits<double>::quiet_NaN()};
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
    if (!line.empty() && line.back() == ',')
        fields.emplace_back();
    return fields;
}

double parseNumber(
    const std::vector<std::string> &fields,
    const std::unordered_map<std::string, std::size_t> &column,
    const std::string &name)
{
    const auto it = column.find(name);
    if (it == column.end() || it->second >= fields.size())
        throw std::runtime_error("Missing benchmark CSV column: " + name);
    return std::stod(fields[it->second]);
}

std::string parseText(
    const std::vector<std::string> &fields,
    const std::unordered_map<std::string, std::size_t> &column,
    const std::string &name)
{
    const auto it = column.find(name);
    if (it == column.end() || it->second >= fields.size())
        throw std::runtime_error("Missing benchmark CSV column: " + name);
    return fields[it->second];
}

std::vector<Request> readRequests(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open benchmark request CSV: " + path.string());

    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("Benchmark request CSV is empty.");
    const auto header = split(line);
    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t i = 0; i < header.size(); ++i)
        column.emplace(header[i], i);

    std::vector<Request> requests;
    while (std::getline(input, line))
    {
        if (line.empty())
            continue;
        const auto fields = split(line);
        Request request;
        request.id = parseText(fields, column, "id");
        request.model = parseText(fields, column, "model");
        request.operation = parseText(fields, column, "operation");
        request.temperature = parseNumber(fields, column, "temperature_K");
        request.pressure = parseNumber(fields, column, "pressure_Pa");
        request.salinity = parseNumber(fields, column, "salinity_molal");
        request.z0 = parseNumber(fields, column, "z0");
        request.z1 = parseNumber(fields, column, "z1");
        requests.push_back(std::move(request));
    }
    return requests;
}

OgEos makePrEos()
{
    // ThermoPack 2.2.3 C1/C2 database values, read through its public API.
    const std::array<double, 2> tc{190.555, 305.4};
    const std::array<double, 2> pc{4.598837e6, 4.883900e6};
    const std::array<double, 2> vc{9.775301080859746e-5, 1.4817717640030296e-4};
    const std::array<double, 2> omega{0.01131, 0.09800};
    const std::array<double, 2> mw{0.016043, 0.030070};
    const std::array<std::array<double, 2>, 2> kij{};
    return OgEos(
        0.4572355289213822, 0.07779607390388846,
        MPMC::CompositionalMixture<OgIndices>(tc, pc, vc, omega, mw, kij),
        1, 2.414213562373095, -0.414213562373095, 1.0e-30);
}

OgEos makeCpaEos()
{
    // Tc/Pc/omega only complete the mixture metadata.  The CPA cubic and
    // association terms below use ThermoPack's explicit H2O/MEOH parameters.
    // These are the component-database values used by ThermoPack's Classic
    // CPA alpha function (not modern reference-property critical constants).
    const std::array<double, 2> tc{647.3, 512.60};
    const std::array<double, 2> pc{22.0483e6, 8.0959e6};
    const std::array<double, 2> vc{5.60e-5, 1.18e-4};
    const std::array<double, 2> omega{0.3443, 0.5650};
    const std::array<double, 2> mw{0.01801528, 0.03204186};
    const std::array<std::array<double, 2>, 2> kij{{
        {{0.0, -0.09}},
        {{-0.09, 0.0}}
    }};
    OgEos eos(
        0.42748, 0.08664,
        MPMC::CompositionalMixture<OgIndices>(tc, pc, vc, omega, mw, kij),
        1, 1.0, 0.0, 1.0e-30);

    OgEos::CubicPlusAssociationOptions cpa;
    cpa.a0 = {0.12277, 0.40531};
    cpa.b = {1.4515e-5, 3.0978e-5};
    cpa.c1 = {0.67359, 0.43102};
    cpa.associationEnergy = {16655.0, 24591.0};
    cpa.associationVolume = {0.06920, 0.01610};
    cpa.donorSites = {2, 1};
    cpa.acceptorSites = {2, 1};
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

SwEos makeSwEos(double salinity, bool useChabab2019 = true)
{
    // Component order is CO2/H2O to make the reported dissolved-CO2 value z0.
    const std::array<double, 2> tc{304.1282, 647.096};
    const std::array<double, 2> pc{7.3773e6, 22.064e6};
    const std::array<double, 2> vc{9.40e-5, 5.60e-5};
    const std::array<double, 2> omega{0.22394, 0.34430};
    const std::array<double, 2> mw{0.0440095, 0.01801528};
    const std::array<std::array<double, 2>, 2> kij{{
        {{0.0, 0.1896}},
        {{0.1896, 0.0}}
    }};
    SwEos eos(
        0.4572355289213822, 0.07779607390388846,
        MPMC::CompositionalMixture<SwIndices>(tc, pc, vc, omega, mw, kij),
        5, 2.414213562373095, -0.414213562373095, 1.0e-30);
    SwEos::SoreideWhitsonOptions sw;
    sw.waterComponent = 1;
    sw.salinityMolality = salinity;
    if (useChabab2019)
        sw.aqueousWaterBip[0] = [](double temperature, double molality) {
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(
                temperature, molality);
        };
    else
        sw.aqueousWaterBip[0] = [](double temperature, double molality) {
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
                temperature, 304.1282, molality);
        };
    sw.aqueousWaterBip[1] = [](double, double) { return 0.0; };
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

template <class Result>
Output fromFlash(const std::string &id, const Result &result)
{
    Output output;
    output.id = id;
    output.converged = result.converged;
    output.phaseCode = result.converged ? result.presence.bits() : 0;
    output.iterations = result.iterations;
    output.beta = result.phaseMoleFraction;
    output.x = result.composition;
    output.zFactor = result.compressibility;
    output.molarDensity = result.molarDensity;
    return output;
}

template <class Eos, class Composition>
Output phaseOutput(
    const Request &request,
    const Eos &eos,
    const Composition &composition,
    MPMC::CompositionalPhase phase)
{
    const auto result = eos.phaseResult(
        request.pressure, request.temperature, composition, phase, false);
    Output output;
    output.id = request.id;
    output.converged = true;
    output.phaseCode = phase == MPMC::CompositionalPhase::Gas ? 2 : 1;
    const std::size_t slot = phase == MPMC::CompositionalPhase::Gas ? 1u :
        phase == MPMC::CompositionalPhase::Water ? 2u : 0u;
    output.x[slot] = composition;
    output.zFactor[slot] = result.compressibility;
    output.molarDensity[slot] = eos.molarDensity(
        request.pressure, request.temperature, composition,
        result.compressibility);
    for (std::size_t i = 0; i < 2; ++i)
        output.lnPhi[i] = std::log(result.fugacityCoefficient[i]);
    return output;
}

template <class Flash, class Composition>
auto restrictedFlash(
    const Flash &flash, double pressure, double temperature,
    const Composition &z)
{
    return flash.flashRestricted(
        pressure, temperature, z,
        MPMC::PhasePresence(
            MPMC::PhasePresence::oilBit | MPMC::PhasePresence::gasBit));
}

template <class Eos, class Flash, class Composition>
std::pair<double, decltype(restrictedFlash(
    std::declval<const Flash &>(), 1.0, 1.0,
    std::declval<const Composition &>()))>
bubblePressure(
    const Eos &eos, const Flash &flash, double temperature,
    double referencePressure,
    const Composition &z)
{
    using Result = decltype(restrictedFlash(flash, 1.0, 1.0, z));

    // A pure-component bubble point has no finite composition split for the
    // flash to expose.  Locate its saturation pressure from equality of the
    // liquid and vapour fugacities of the same production EOS.  Above the
    // component critical temperature no bracket exists and NaN is returned.
    const bool pureEndpoint = z[0] <= 1.0e-14 || z[1] <= 1.0e-14;
    if (pureEndpoint)
    {
        const std::size_t component = z[0] > z[1] ? 0u : 1u;
        auto residual = [&](double pressure) {
            const auto liquid = eos.phaseResult(
                pressure, temperature, z, MPMC::CompositionalPhase::Oil, false);
            const auto vapour = eos.phaseResult(
                pressure, temperature, z, MPMC::CompositionalPhase::Gas, false);
            if (std::abs(liquid.compressibility - vapour.compressibility) < 1.0e-7)
                return std::numeric_limits<double>::quiet_NaN();
            return std::log(liquid.fugacityCoefficient[component]) -
                std::log(vapour.fugacityCoefficient[component]);
        };

        constexpr double pMin = 1.0e3;
        constexpr double pMax = 2.0e8;
        double lower = pMin;
        double fLower = residual(lower);
        for (int sample = 1; sample <= 240; ++sample)
        {
            const double fraction = static_cast<double>(sample) / 240.0;
            const double upper = std::exp(
                std::log(pMin) + fraction * (std::log(pMax) - std::log(pMin)));
            const double fUpper = residual(upper);
            if (std::isfinite(fLower) && std::isfinite(fUpper) &&
                fLower * fUpper <= 0.0)
            {
                double a = lower;
                double b = upper;
                double fa = fLower;
                for (int iteration = 0; iteration < 80; ++iteration)
                {
                    const double midpoint = std::sqrt(a * b);
                    const double fm = residual(midpoint);
                    if (!std::isfinite(fm))
                        break;
                    if (fa * fm <= 0.0)
                    {
                        b = midpoint;
                    }
                    else
                    {
                        a = midpoint;
                        fa = fm;
                    }
                    if (b / a - 1.0 < 1.0e-10)
                        break;
                }
                return {std::sqrt(a * b), Result{}};
            }
            lower = upper;
            fLower = fUpper;
        }
        return {std::numeric_limits<double>::quiet_NaN(), Result{}};
    }

    auto findBracket = [&](double pMin, double pMax, int samples)
        -> std::pair<double, double> {
        double previousP = pMin;
        int previousCode = 0;
        for (int i = 0; i < samples; ++i)
        {
            const double f = static_cast<double>(i) / static_cast<double>(samples - 1);
            const double pressure = std::exp(
                std::log(pMin) + f * (std::log(pMax) - std::log(pMin)));
            const auto result = restrictedFlash(flash, pressure, temperature, z);
            const int code = result.converged ? result.presence.bits() : 0;
            if (previousCode == 3 && code == 1)
                return {previousP, pressure};
            if (code != 0)
            {
                previousCode = code;
                previousP = pressure;
            }
        }
        return {0.0, 0.0};
    };

    // Prefer the transition nearest the supplied experimental/continuation
    // pressure.  A very broad first scan can encounter a small numerical
    // two-phase island near the critical region before the physical bubble
    // branch and make a dense P-x curve jump between roots.
    auto bracket = findBracket(
        std::max(1.0e3, referencePressure / 1.35),
        std::min(2.0e8, referencePressure * 1.35), 65);
    if (!(bracket.first > 0.0))
        bracket = findBracket(
            std::max(1.0e3, referencePressure / 4.0),
            std::min(2.0e8, referencePressure * 4.0), 65);
    if (!(bracket.first > 0.0))
        bracket = findBracket(1.0e3, 2.0e8, 129);
    if (!(bracket.first > 0.0))
        return {std::numeric_limits<double>::quiet_NaN(), Result{}};

    Result twoPhase = restrictedFlash(
        flash, bracket.first, temperature, z);
    for (int iteration = 0; iteration < 60; ++iteration)
    {
        const double pressure = std::sqrt(bracket.first * bracket.second);
        const auto result = restrictedFlash(flash, pressure, temperature, z);
        const int code = result.converged ? result.presence.bits() : 0;
        if (code == 3)
        {
            bracket.first = pressure;
            twoPhase = result;
        }
        else if (code == 1)
        {
            bracket.second = pressure;
        }
        else
        {
            break;
        }
        if (bracket.second / bracket.first - 1.0 < 1.0e-9)
            break;
    }
    return {std::sqrt(bracket.first * bracket.second), twoPhase};
}

template <class Eos, class Flash, class Composition>
std::pair<double, decltype(restrictedFlash(
    std::declval<const Flash &>(), 1.0, 1.0,
    std::declval<const Composition &>()))>
bubbleTemperature(
    const Eos &eos, const Flash &flash, double pressure,
    double referenceTemperature,
    const Composition &z)
{
    using Result = decltype(restrictedFlash(flash, 1.0, 1.0, z));

    // At a pure-component endpoint the two-phase flash has no finite-width
    // composition split to expose, although the saturation temperature is
    // still defined.  Locate it directly from equality of pure liquid/vapour
    // fugacities.  This is the same EOS criterion used by the flash, not an
    // experimental-data fit or endpoint interpolation.
    const bool pureEndpoint = z[0] <= 1.0e-14 || z[1] <= 1.0e-14;
    if (pureEndpoint)
    {
        const std::size_t component = z[0] > z[1] ? 0u : 1u;
        auto residual = [&](double temperature) {
            const auto liquid = eos.phaseResult(
                pressure, temperature, z, MPMC::CompositionalPhase::Oil, false);
            const auto vapour = eos.phaseResult(
                pressure, temperature, z, MPMC::CompositionalPhase::Gas, false);
            if (std::abs(liquid.compressibility - vapour.compressibility) < 1.0e-7)
                return std::numeric_limits<double>::quiet_NaN();
            return std::log(liquid.fugacityCoefficient[component]) -
                std::log(vapour.fugacityCoefficient[component]);
        };
        double lower = 250.0;
        double fLower = residual(lower);
        for (int sample = 1; sample <= 160; ++sample)
        {
            const double upper = 250.0 + 400.0 * sample / 160.0;
            const double fUpper = residual(upper);
            if (std::isfinite(fLower) && std::isfinite(fUpper) &&
                fLower * fUpper <= 0.0)
            {
                double a = lower;
                double b = upper;
                double fa = fLower;
                for (int iteration = 0; iteration < 80; ++iteration)
                {
                    const double midpoint = 0.5 * (a + b);
                    const double fm = residual(midpoint);
                    if (fa * fm <= 0.0)
                    {
                        b = midpoint;
                    }
                    else
                    {
                        a = midpoint;
                        fa = fm;
                    }
                    if (b - a < 1.0e-8)
                        break;
                }
                return {0.5 * (a + b), Result{}};
            }
            lower = upper;
            fLower = fUpper;
        }
        return {std::numeric_limits<double>::quiet_NaN(), Result{}};
    }

    auto findBracket = [&](double tMin, double tMax, int samples)
        -> std::pair<double, double> {
        double previousT = tMin;
        int previousCode = 0;
        for (int i = 0; i < samples; ++i)
        {
            const double f = static_cast<double>(i) / static_cast<double>(samples - 1);
            const double temperature = tMin + f * (tMax - tMin);
            const auto result = restrictedFlash(flash, pressure, temperature, z);
            const int code = result.converged ? result.presence.bits() : 0;
            if (previousCode == 1 && code == 3)
                return {previousT, temperature};
            if (code != 0)
            {
                previousCode = code;
                previousT = temperature;
            }
        }
        return {0.0, 0.0};
    };

    auto bracket = findBracket(
        std::max(250.0, referenceTemperature - 35.0),
        std::min(650.0, referenceTemperature + 35.0), 29);
    if (!(bracket.first > 0.0))
        bracket = findBracket(250.0, 650.0, 81);
    if (!(bracket.first > 0.0))
        return {std::numeric_limits<double>::quiet_NaN(), Result{}};

    Result twoPhase = restrictedFlash(
        flash, pressure, bracket.second, z);
    for (int iteration = 0; iteration < 60; ++iteration)
    {
        const double temperature = 0.5 * (bracket.first + bracket.second);
        const auto result = restrictedFlash(flash, pressure, temperature, z);
        const int code = result.converged ? result.presence.bits() : 0;
        if (code == 1)
        {
            bracket.first = temperature;
        }
        else if (code == 3)
        {
            bracket.second = temperature;
            twoPhase = result;
        }
        else
        {
            break;
        }
        if (bracket.second - bracket.first < 1.0e-7)
            break;
    }
    return {0.5 * (bracket.first + bracket.second), twoPhase};
}

Output runOgRequest(const Request &request, const OgEos &eos)
{
    const OgComposition composition{request.z0, request.z1};
    if (request.operation == "phase_liquid")
        return phaseOutput(
            request, eos, composition, MPMC::CompositionalPhase::Oil);
    if (request.operation == "phase_vapor")
        return phaseOutput(
            request, eos, composition, MPMC::CompositionalPhase::Gas);

    MPMC::ThreePhaseFlashOptions options;
    options.compositionFloor = 1.0e-30;
    const OgFlash flash(eos, options);
    if (request.operation == "flash")
        return fromFlash(
            request.id,
            restrictedFlash(
                flash, request.pressure, request.temperature, composition));
    if (request.operation == "bubble_pressure")
    {
        auto [value, result] = bubblePressure(
            eos, flash, request.temperature, request.pressure, composition);
        if ((composition[0] <= 1.0e-14 || composition[1] <= 1.0e-14) &&
            std::isfinite(value))
        {
            Output output;
            output.id = request.id;
            output.converged = true;
            output.phaseCode = 3;
            output.x[0] = composition;
            output.x[1] = composition;
            output.value = value;
            return output;
        }
        Output output = fromFlash(request.id, result);
        output.value = value;
        return output;
    }
    if (request.operation == "bubble_temperature")
    {
        auto [value, result] = bubbleTemperature(
            eos, flash, request.pressure, request.temperature, composition);
        Output output = fromFlash(request.id, result);
        output.value = value;
        return output;
    }
    throw std::runtime_error("Unsupported O/G benchmark operation: " + request.operation);
}

Output runSwRequest(const Request &request)
{
    if (request.operation == "sw_water_alpha")
    {
        Output output;
        output.id = request.id;
        output.converged = true;
        output.value = MPMC::SoreideWhitsonCorrelations::waterAlpha(
            request.temperature, 647.096, request.salinity);
        return output;
    }
    if (request.operation == "sw_co2_bip")
    {
        Output output;
        output.id = request.id;
        output.converged = true;
        output.value = MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
            request.temperature, 304.1282, request.salinity);
        return output;
    }
    if (request.operation == "sw_co2_bip_chabab2019")
    {
        Output output;
        output.id = request.id;
        output.converged = true;
        output.value = MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(
            request.temperature, request.salinity);
        return output;
    }

    const bool useChabab2019 = request.operation != "sw_solubility_legacy";
    const auto eos = makeSwEos(request.salinity, useChabab2019);
    const SwComposition composition{request.z0, request.z1};
    if (request.operation == "phase_water")
        return phaseOutput(
            request, eos, composition, MPMC::CompositionalPhase::Water);
    if (request.operation == "phase_vapor")
        return phaseOutput(
            request, eos, composition, MPMC::CompositionalPhase::Gas);

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 1;
    const SwFlash flash(eos, options);
    const auto result = flash.flash(
        request.pressure, request.temperature, composition);
    Output output = fromFlash(request.id, result);
    if ((request.operation == "sw_solubility" ||
         request.operation == "sw_solubility_legacy") && result.converged)
    {
        const double xCo2 = result.composition[2][0];
        const double xWater = result.composition[2][1];
        output.value = xCo2 / (xWater * 0.01801528);
    }
    else if (request.operation != "flash" &&
             request.operation != "sw_solubility" &&
             request.operation != "sw_solubility_legacy")
    {
        throw std::runtime_error("Unsupported SW benchmark operation: " + request.operation);
    }
    return output;
}

void writeOutputs(
    const std::filesystem::path &path,
    const std::vector<Output> &outputs)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot open benchmark output CSV: " + path.string());
    out << "id,converged,phase_code,iterations,beta_o,beta_g,beta_w,"
           "xo0,xo1,yg0,yg1,xw0,xw1,Z_o,Z_g,Z_w,rho_o,rho_g,rho_w,"
           "lnphi0,lnphi1,value\n";
    out << std::setprecision(17);
    for (const auto &row : outputs)
    {
        out << row.id << ',' << (row.converged ? 1 : 0) << ','
            << row.phaseCode << ',' << row.iterations;
        for (double value : row.beta)
            out << ',' << value;
        for (const auto &phase : row.x)
            for (double value : phase)
                out << ',' << value;
        for (double value : row.zFactor)
            out << ',' << value;
        for (double value : row.molarDensity)
            out << ',' << value;
        out << ',' << row.lnPhi[0] << ',' << row.lnPhi[1]
            << ',' << row.value << '\n';
    }
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            std::cerr << "Usage: three_eos_public_benchmark requests.csv results.csv\n";
            return 2;
        }
        const auto requests = readRequests(argv[1]);
        std::vector<Output> outputs;
        outputs.reserve(requests.size());
        const OgEos pr = makePrEos();
        const OgEos cpa = makeCpaEos();
        for (const auto &request : requests)
        {
            try
            {
                if (request.model == "pr")
                    outputs.push_back(runOgRequest(request, pr));
                else if (request.model == "cpa")
                    outputs.push_back(runOgRequest(request, cpa));
                else if (request.model == "sw")
                    outputs.push_back(runSwRequest(request));
                else
                    throw std::runtime_error("Unsupported benchmark model: " + request.model);
            }
            catch (const std::exception &error)
            {
                Output failed;
                failed.id = request.id;
                outputs.push_back(failed);
                std::cerr << "[three-eos-public-benchmark] row " << request.id
                          << " failed: " << error.what() << '\n';
            }
        }
        writeOutputs(argv[2], outputs);
        std::cout << "[three-eos-public-benchmark] wrote " << outputs.size()
                  << " rows to " << argv[2] << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[three-eos-public-benchmark] ERROR: "
                  << error.what() << '\n';
        return 1;
    }
}
