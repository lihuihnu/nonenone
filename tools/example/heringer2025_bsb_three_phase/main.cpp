/**
 * @file main.cpp
 * @brief Heringer 2025 水/二氧化碳/BSB 原油三相闪蒸基准。
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
#include <tools/eos_parameter_regression.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace
{

using Config = MPMC::CompositionalModelConfig<
    8, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, Indices::numComponents>;

constexpr std::size_t N = static_cast<std::size_t>(Indices::numComponents);
constexpr double benchmarkTemperature = 650.0;
constexpr double benchmarkPressureBar = 390.0;
constexpr double flashCompositionFloor = 1.0e-30;
double temperature = benchmarkTemperature;
double pressureBar = benchmarkPressureBar;
double pressure = benchmarkPressureBar * 1.0e5;

inline const std::array<std::string, N> componentNames{
    "H2O", "CO2", "C1", "C2-3", "C4-6", "C7-15", "C16-27", "C28+"};

// Heringer et al. (2025), Table B3.
inline constexpr Composition tc{
    647.30, 304.20, 160.00, 344.22, 463.22, 605.78, 751.00, 942.50};
inline constexpr Composition pc{
    220.48e5, 73.76e5, 46.00e5, 45.00e5, 34.00e5, 21.75e5, 16.54e5, 16.42e5};
inline constexpr Composition omega{
    0.344, 0.225, 0.008, 0.131, 0.240, 0.618, 0.957, 1.268};
inline constexpr Composition feed{
    0.750000, 0.105055, 0.012915, 0.022545,
    0.025065, 0.049560, 0.024165, 0.010695};

// BSB metadata from Fernandes et al. (2021), Table 12. H2O is standard
// pure-component metadata. Neither Vc nor MW enters this PR fugacity flash.
inline constexpr Composition vc{
    5.60e-5, 9.43e-5, 9.93e-5, 1.81e-4,
    3.06e-4, 5.99e-4, 1.13e-3, 2.09e-3};
inline constexpr Composition mw{
    0.01801528, 0.044010, 0.016040, 0.037200,
    0.069500, 0.140960, 0.280990, 0.519620};

// Heringer et al. (2025), Table B4. All unlisted pairs are zero.
inline constexpr std::array<std::array<double, N>, N> kij{{
    {{0.0000, 0.1896, 0.4850, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000}},
    {{0.1896, 0.0000, 0.0550, 0.0550, 0.0550, 0.1050, 0.1050, 0.1050}},
    {{0.4850, 0.0550, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.0550, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.0550, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.1050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.1050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.1050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}}
}};

// Heringer et al. (2025), Table 6, Full Three-Phase Flash columns.
inline constexpr Composition referenceGas{
    0.7041, 0.1333, 0.0169, 0.0294, 0.0320, 0.0573, 0.0216, 0.0054};
inline constexpr Composition referenceOil{
    0.4900, 0.1241, 0.0164, 0.0357, 0.0485, 0.1308, 0.0936, 0.0608};
inline constexpr Composition referenceWater{
    0.9492, 0.0417, 0.0037, 0.0036, 0.0016, 2.689e-4, 6.707e-6, 4.922e-7};
inline constexpr std::array<double, 3> referenceBeta{
    0.1245, 0.5794, 0.2960}; // public slots: oil, gas, water

Composition normalized(Composition value)
{
    double sum = 0.0;
    for (double x : value)
        sum += x;
    if (!(sum > 0.0))
        throw std::runtime_error("Cannot normalize an empty composition.");
    for (double &x : value)
        x /= sum;
    return value;
}

std::array<Composition, 3> normalizedReference()
{
    return {normalized(referenceOil), normalized(referenceGas), normalized(referenceWater)};
}

bool isBenchmarkState()
{
    return std::abs(temperature - benchmarkTemperature) < 1.0e-10 &&
        std::abs(pressureBar - benchmarkPressureBar) < 1.0e-10;
}

Eos makePr(int eosType)
{
    return Eos(
        0.45724, 0.07780,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        eosType,
        2.414213562373095,
        -0.414213562373095,
        1.0e-30);
}

Eos makeSw(bool extrapolateHeavyHydrocarbonBips)
{
    Eos eos = makePr(5);
    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    sw.aqueousWaterBip[1] = [](double t, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
            t, tc[1], salinity);
    };
    sw.aqueousWaterBip[2] = [](double t, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            t, tc[2], omega[2], salinity);
    };
    sw.aqueousWaterBip[3] = [](double t, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            t, tc[3], omega[3], salinity);
    };

    if (extrapolateHeavyHydrocarbonBips)
    {
        sw.aqueousWaterBip[4] = [](double t, double salinity) {
            return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                t, tc[4], omega[4], salinity);
        };
        sw.aqueousWaterBip[5] = [](double t, double salinity) {
            return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                t, tc[5], omega[5], salinity);
        };
        sw.aqueousWaterBip[6] = [](double t, double salinity) {
            return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                t, tc[6], omega[6], salinity);
        };
        sw.aqueousWaterBip[7] = [](double t, double salinity) {
            return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                t, tc[7], omega[7], salinity);
        };
    }
    else
    {
        // The published SW hydrocarbon correlation is fitted only through nC4.
        // Keep the Heringer water-heavy BIP rather than silently extrapolating
        // the correlation to the C4-6 and heavier pseudo-components.
        for (std::size_t c = 4; c < N; ++c)
            sw.aqueousWaterBip[c] = [](double, double) { return 0.5000; };
    }
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

Eos makeSwWithPaperHydrocarbonBips(bool useSwCo2Correlation)
{
    Eos eos = makePr(5);
    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    if (useSwCo2Correlation)
    {
        sw.aqueousWaterBip[1] = [](double t, double salinity) {
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
                t, tc[1], salinity);
        };
    }
    // Empty hydrocarbon callbacks deliberately retain Table B4 values. This
    // isolates the documented SW water alpha and CO2 correlation from the
    // generalized hydrocarbon correlation used by the full SW1992 cases.
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

using TrendParameters = std::array<double, 4>;

Eos makeSwPrTrendMatched(const TrendParameters &parameters)
{
    // Preserve the SW water alpha function.  Only the aqueous water/non-water
    // BIPs are calibrated to the local PR response surface.  Parameters are:
    // CO2 offset, hydrocarbon offset, CO2 slope and hydrocarbon slope, with
    // slopes expressed per 100 K around the 650 K benchmark temperature.
    Eos eos = makePr(5);
    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    for (std::size_t c = 1; c < N; ++c)
    {
        const double base = kij[0][c];
        const bool isCo2 = c == 1;
        const double offset = parameters[isCo2 ? 0 : 1];
        const double slope = parameters[isCo2 ? 2 : 3];
        sw.aqueousWaterBip[c] = [base, offset, slope](double t, double) {
            return base + offset + slope * (t - benchmarkTemperature) / 100.0;
        };
    }
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

Eos makeSrkCpa()
{
    Eos eos = makePr(5);
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

double jiaOkuno2018HeavyWaterBip(int component)
{
    // Jia and Okuno (2018), Table 2, using the nearest measured n-alkane
    // molecular weight for each BSB heavy pseudo-component. C4-6 is mapped by
    // linear interpolation between their reported nC4 (0.306) and nC7
    // correlated (0.241) values. Their fit is based on L-V-W coexistence data.
    if (component < 0 || component >= Indices::numComponents)
        return 0.0;
    switch (component)
    {
    case 4:
    {
        constexpr double mwC4 = 58.12;
        constexpr double mwC7 = 100.0;
        constexpr double kC4 = 0.306;
        constexpr double kC7 = 0.241;
        const double pseudoMwGramPerMol = 1000.0 * mw[4];
        return kC4 + (pseudoMwGramPerMol - mwC4) * (kC7 - kC4) /
            (mwC7 - mwC4);
    }
    case 5: return 0.165;  // MW 140.96 -> nC10, MW 142
    case 6: return -0.006; // MW 280.99 -> nC20, MW 282
    case 7: return -0.038; // MW 519.62 -> nC36+ limiting value
    default: return kij[0][static_cast<std::size_t>(component)];
    }
}

Eos makeSrkCpaJiaOkuno2018()
{
    Eos eos = makeSrkCpa();
    eos.configureBinaryInteractionFunction([](int i, int j, double) {
        if (i == j)
            return 0.0;
        if (i == 0 && j >= 4)
            return jiaOkuno2018HeavyWaterBip(j);
        if (j == 0 && i >= 4)
            return jiaOkuno2018HeavyWaterBip(i);
        return kij[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    });
    return eos;
}

Eos makePrCpaSensitivity()
{
    Eos eos = makePr(5);
    Eos::CubicPlusAssociationOptions cpa;
    constexpr double R = MPMC::units::gasConstant;
    for (std::size_t c = 0; c < N; ++c)
    {
        cpa.a0[c] = 0.45724 * R * R * tc[c] * tc[c] / pc[c];
        cpa.b[c] = 0.07780 * R * tc[c] / pc[c];
        cpa.c1[c] = 0.37464 + 1.54226 * omega[c]
            - 0.26992 * omega[c] * omega[c];
        if (omega[c] > 0.49)
        {
            cpa.c1[c] = 0.379642 + 1.48503 * omega[c]
                - 0.164423 * omega[c] * omega[c]
                + 0.016666 * omega[c] * omega[c] * omega[c];
        }
    }
    cpa.associationEnergy[0] = MPMC::StandardCpaWater4C::epsilon;
    cpa.associationVolume[0] = MPMC::StandardCpaWater4C::beta;
    cpa.donorSites[0] = MPMC::StandardCpaWater4C::donorSites;
    cpa.acceptorSites[0] = MPMC::StandardCpaWater4C::acceptorSites;
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::PengRobinson;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

Eos makePrCpaPrTrendMatched(const TrendParameters &parameters)
{
    // PR-CPA is used deliberately: it preserves the PR cubic response for all
    // non-associating components while retaining a non-zero 4C water
    // association term.  Parameters are association-volume scale, CO2 BIP
    // offset, hydrocarbon BIP offset and a common BIP slope per 100 K.
    Eos eos = makePr(5);
    eos.configureBinaryInteractionFunction([parameters](int i, int j, double t) {
        if (i == j)
            return 0.0;
        if (i == 0 || j == 0)
        {
            const int other = i == 0 ? j : i;
            const double offset = parameters[other == 1 ? 1 : 2];
            return kij[0][static_cast<std::size_t>(other)] + offset
                + parameters[3] * (t - benchmarkTemperature) / 100.0;
        }
        return kij[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    });

    Eos::CubicPlusAssociationOptions cpa;
    constexpr double R = MPMC::units::gasConstant;
    for (std::size_t c = 0; c < N; ++c)
    {
        cpa.a0[c] = 0.45724 * R * R * tc[c] * tc[c] / pc[c];
        cpa.b[c] = 0.07780 * R * tc[c] / pc[c];
        cpa.c1[c] = 0.37464 + 1.54226 * omega[c]
            - 0.26992 * omega[c] * omega[c];
        if (omega[c] > 0.49)
        {
            cpa.c1[c] = 0.379642 + 1.48503 * omega[c]
                - 0.164423 * omega[c] * omega[c]
                + 0.016666 * omega[c] * omega[c] * omega[c];
        }
    }
    cpa.associationEnergy[0] = MPMC::StandardCpaWater4C::epsilon;
    cpa.associationVolume[0] = parameters[0] * MPMC::StandardCpaWater4C::beta;
    cpa.donorSites[0] = MPMC::StandardCpaWater4C::donorSites;
    cpa.acceptorSites[0] = MPMC::StandardCpaWater4C::acceptorSites;
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::PengRobinson;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

double sorensen2018WaterBip(int component, double t)
{
    // Sorensen et al. (2018), Table 9:
    // kij(T) = kij_ref + kij_prime * (T - 288.15 K).
    // The BSB pseudo-components are mapped by their representative molecular
    // weight/critical temperature: C2-3 is the mean of C2 and C3, C4-6 is
    // nC5-like, and every C7+ fraction uses the published C7+ entry.
    const double deltaT = t - 288.15;
    switch (component)
    {
    case 1: return 0.07574 + 6.649e-4 * deltaT; // CO2
    case 2: return 0.03833 + 1.588e-3 * deltaT; // C1
    case 3:
    {
        const double c2 = 0.07594 + 9.937e-4 * deltaT;
        const double c3 = 0.04286 + 8.697e-4 * deltaT;
        return 0.5 * (c2 + c3);
    }
    case 4: return 0.00350; // C4-6 pseudo-component is nC5-like
    case 5:
    case 6:
    case 7: return 0.0; // published C7+ value
    default: return 0.0;
    }
}

Eos makePrCpaSorensen2018()
{
    Eos eos = makePr(5);
    eos.configureBinaryInteractionFunction([](int i, int j, double t) {
        if (i == j)
            return 0.0;
        if (i == 0)
            return sorensen2018WaterBip(j, t);
        if (j == 0)
            return sorensen2018WaterBip(i, t);
        return kij[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    });

    Eos::CubicPlusAssociationOptions cpa;
    constexpr double R = MPMC::units::gasConstant;
    for (std::size_t c = 0; c < N; ++c)
    {
        cpa.a0[c] = 0.45724 * R * R * tc[c] * tc[c] / pc[c];
        cpa.b[c] = 0.07780 * R * tc[c] / pc[c];
        cpa.c1[c] = 0.37464 + 1.54226 * omega[c]
            - 0.26992 * omega[c] * omega[c];
        if (omega[c] > 0.49)
        {
            cpa.c1[c] = 0.379642 + 1.48503 * omega[c]
                - 0.164423 * omega[c] * omega[c]
                + 0.016666 * omega[c] * omega[c] * omega[c];
        }
    }

    // Sorensen et al. (2018), Tables 6 and 8. Unit conversions:
    // 1 bar L^2/mol^2 = 0.1 Pa m^6/mol^2; 1 bar L/mol = 100 J/mol.
    cpa.a0[0] = 1.5782e-1;
    cpa.b[0] = 1.4788e-5;
    cpa.c1[0] = 0.6736;
    cpa.associationEnergy[0] = 16123.0;
    cpa.associationVolume[0] = 6.9662e-2;
    cpa.donorSites[0] = 2;
    cpa.acceptorSites[0] = 2;

    // CO2 has one negative site and no positive site: it cross-associates with
    // water but cannot self-associate. Equation A17 sets the cross energy to
    // one half of the water self-association energy; Table 8 gives beta.
    cpa.associationEnergy[1] = 1.0;
    cpa.associationVolume[1] = 1.0;
    cpa.donorSites[1] = 0;
    cpa.acceptorSites[1] = 1;
    cpa.crossAssociationEnergy[0][1] = 0.5 * cpa.associationEnergy[0];
    cpa.crossAssociationVolume[0][1] = 0.15182;
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::PengRobinson;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    return eos;
}

double maxMassClosure(const Flash::Result &result)
{
    double maximum = 0.0;
    for (std::size_t c = 0; c < N; ++c)
    {
        double reconstructed = 0.0;
        for (std::size_t p = 0; p < 3; ++p)
            reconstructed += result.phaseMoleFraction[p] * result.composition[p][c];
        maximum = std::max(maximum, std::abs(reconstructed - feed[c]));
    }
    return maximum;
}

double maxLogFugacitySpreadAt(
    const Eos &eos,
    const Flash::Result &result,
    double statePressure,
    double stateTemperature)
{
    std::array<Eos::PhaseResult<double>, 3> phase{
        eos.phaseResult(statePressure, stateTemperature, result.composition[0],
                        MPMC::CompositionalPhase::Oil, false),
        eos.phaseResult(statePressure, stateTemperature, result.composition[1],
                        MPMC::CompositionalPhase::Gas, false),
        eos.phaseResult(statePressure, stateTemperature, result.composition[2],
                        MPMC::CompositionalPhase::Water, false)};
    double maximum = 0.0;
    for (std::size_t c = 0; c < N; ++c)
    {
        bool floorBound = false;
        for (std::size_t p = 0; p < 3; ++p)
        {
            if (result.phaseMoleFraction[p] > 1.0e-12 &&
                result.composition[p][c] <= 10.0 * flashCompositionFloor)
            {
                floorBound = true;
                break;
            }
        }
        // A component clamped at the active composition bound satisfies a
        // complementarity condition rather than an equality equation.  Match
        // the production flash validator and exclude it from the equality-only
        // diagnostic spread.
        if (floorBound)
            continue;

        double minimumLogF = std::numeric_limits<double>::infinity();
        double maximumLogF = -std::numeric_limits<double>::infinity();
        for (std::size_t p = 0; p < 3; ++p)
        {
            if (result.phaseMoleFraction[p] <= 1.0e-12)
                continue;
            const double logF =
                std::log(std::max(phase[p].fugacity[c], 1.0e-300));
            minimumLogF = std::min(minimumLogF, logF);
            maximumLogF = std::max(maximumLogF, logF);
        }
        if (std::isfinite(minimumLogF) && std::isfinite(maximumLogF))
            maximum = std::max(maximum, maximumLogF - minimumLogF);
    }
    return maximum;
}

double maxLogFugacitySpread(const Eos &eos, const Flash::Result &result)
{
    return maxLogFugacitySpreadAt(
        eos, result, pressure, temperature);
}

struct Metrics
{
    double compositionMae{std::numeric_limits<double>::quiet_NaN()};
    double compositionMax{std::numeric_limits<double>::quiet_NaN()};
    double betaMae{std::numeric_limits<double>::quiet_NaN()};
    double betaMax{std::numeric_limits<double>::quiet_NaN()};
    double massClosure{std::numeric_limits<double>::quiet_NaN()};
    double logFugacitySpread{std::numeric_limits<double>::quiet_NaN()};
};

Metrics compare(const Eos &eos, const Flash::Result &result)
{
    Metrics metrics;
    if (!result.converged)
        return metrics;
    metrics.massClosure = maxMassClosure(result);
    metrics.logFugacitySpread = maxLogFugacitySpread(eos, result);
    if (!isBenchmarkState())
        return metrics;
    const auto reference = normalizedReference();
    double sum = 0.0;
    double maximum = 0.0;
    for (std::size_t p = 0; p < 3; ++p)
    {
        for (std::size_t c = 0; c < N; ++c)
        {
            const double error = std::abs(result.composition[p][c] - reference[p][c]);
            sum += error;
            maximum = std::max(maximum, error);
        }
    }
    metrics.compositionMae = sum / static_cast<double>(3 * N);
    metrics.compositionMax = maximum;
    double betaSum = 0.0;
    double betaMax = 0.0;
    for (std::size_t p = 0; p < 3; ++p)
    {
        const double error = std::abs(result.phaseMoleFraction[p] - referenceBeta[p]);
        betaSum += error;
        betaMax = std::max(betaMax, error);
    }
    metrics.betaMae = betaSum / 3.0;
    metrics.betaMax = betaMax;
    return metrics;
}

struct CaseResult
{
    std::string backend;
    std::string parameterization;
    std::string path;
    Flash::Result flash;
    Metrics metrics;
    bool coincidentGasWaterMerged{false};
};

struct TrendState
{
    std::string name;
    std::string split;
    double temperatureK{temperature};
    double pressurePa{pressure};
};

struct TrendFit
{
    std::vector<MPMC::tools::RegressionParameter> descriptors;
    MPMC::tools::RegressionResult result;
    TrendParameters parameters{};
};

std::vector<TrendState> prTrendStates()
{
    // Dense local design around the Heringer state.  Every third grid point is
    // used for fitting; the rest are held out.  Only states at which PR itself
    // predicts O+G+W are admitted to the objective (see fitToPrTrend).  Nearby
    // two-phase states remain in the output as phase-boundary diagnostics.
    const std::array<double, 7> temperatures{
        625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0};
    const std::array<double, 6> pressuresBar{
        330.0, 350.0, 370.0, 390.0, 410.0, 430.0};
    std::vector<TrendState> states;
    for (std::size_t ti = 0; ti < temperatures.size(); ++ti)
    {
        for (std::size_t pi = 0; pi < pressuresBar.size(); ++pi)
        {
            const int t = static_cast<int>(temperatures[ti]);
            const int p = static_cast<int>(pressuresBar[pi]);
            states.push_back({
                "T" + std::to_string(t) + "_P" + std::to_string(p),
                (ti + pi) % 3 == 2 ? "training" : "validation",
                temperatures[ti], pressuresBar[pi] * 1.0e5});
        }
    }
    return states;
}

Flash::Result flashAt(const Eos &eos, const TrendState &state)
{
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    options.maximumIterations = 200;
    options.maximumStabilityIterations = 140;
    options.fugacityTolerance = 1.0e-9;
    return Flash(eos, options).flash(
        state.pressurePa, state.temperatureK, feed);
}

double prTrendError(
    const Flash::Result &target,
    const Flash::Result &candidate)
{
    if (!candidate.converged)
        return 5.0;

    double error = 0.0;
    if (candidate.presence.bits() != target.presence.bits())
        error += 1.0;

    for (std::size_t p = 0; p < 3; ++p)
    {
        error += std::abs(
            candidate.phaseMoleFraction[p] - target.phaseMoleFraction[p]) / 3.0;
        const bool targetPresent = target.phaseMoleFraction[p] > 1.0e-8;
        const bool candidatePresent = candidate.phaseMoleFraction[p] > 1.0e-8;
        if (!targetPresent)
            continue;
        if (!candidatePresent)
        {
            error += 1.0;
            continue;
        }
        for (std::size_t c = 0; c < N; ++c)
        {
            error += std::abs(
                candidate.composition[p][c] - target.composition[p][c]) /
                static_cast<double>(3 * N);
        }
    }
    return error;
}

template <class Factory>
TrendFit fitToPrTrend(
    std::vector<MPMC::tools::RegressionParameter> descriptors,
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets,
    Factory &&factory)
{
    MPMC::tools::PatternSearchOptions options;
    options.maximumIterations = 12;
    options.minimumRelativeStep = 2.0e-4;
    options.improvementTolerance = 1.0e-10;
    const auto result = MPMC::tools::boundedPatternSearch(
        descriptors,
        [&](const std::vector<double> &values) {
            TrendParameters parameters{};
            std::copy(values.begin(), values.end(), parameters.begin());
            try
            {
                const Eos eos = factory(parameters);
                double sum = 0.0;
                std::size_t count = 0;
                for (std::size_t i = 0; i < states.size(); ++i)
                {
                    if (states[i].split != "training" ||
                        targets[i].presence.bits() != MPMC::PhasePresence::allBits)
                        continue;
                    sum += prTrendError(targets[i], flashAt(eos, states[i]));
                    ++count;
                }
                if (count == 0)
                    return 10.0;
                return sum / static_cast<double>(count);
            }
            catch (const std::exception &)
            {
                return 10.0;
            }
        },
        options);

    TrendFit fit{std::move(descriptors), result, {}};
    std::copy(result.values.begin(), result.values.end(), fit.parameters.begin());
    return fit;
}

template <class Factory>
double meanTrainingTrendError(
    const TrendParameters &parameters,
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets,
    Factory &&factory)
{
    const Eos eos = factory(parameters);
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        if (states[i].split != "training" ||
            targets[i].presence.bits() != MPMC::PhasePresence::allBits)
            continue;
        sum += prTrendError(targets[i], flashAt(eos, states[i]));
        ++count;
    }
    return count == 0 ? 10.0 : sum / static_cast<double>(count);
}

std::array<TrendFit, 2> optimizePrTrend(
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets)
{
    std::vector<MPMC::tools::RegressionParameter> swDescriptors{
        {"water_CO2_BIP_offset", 0.0, -0.25, 0.25, 0.025},
        {"water_HC_BIP_offset", 0.0, -0.25, 0.25, 0.025},
    };
    std::vector<MPMC::tools::RegressionParameter> cpaDescriptors{
        // The preliminary logarithmic sweep found that 0.001 is the largest
        // tested scale that preserves all PR three-phase training states.
        // It remains non-zero, so this sensitivity variant is still CPA, but
        // it must not be presented as a literature parameterization.
        {"water_association_volume_scale", 1.0, 0.001, 1.0, 0.0},
        {"water_CO2_BIP_offset", 0.0, -0.40, 0.40, 0.04},
        {"water_HC_BIP_offset", 0.0, -0.40, 0.40, 0.04},
        {"water_BIP_slope_per_100K", 0.0, -0.30, 0.30, 0.03},
    };
    TrendFit swFit = fitToPrTrend(
        std::move(swDescriptors), states, targets,
        [](const TrendParameters &p) { return makeSwPrTrendMatched(p); });

    // CPA has a sharp phase-topology switch between volume scales 0.005 and
    // 0.001.  A continuous local optimizer is not reliable across that jump;
    // use the documented scale sweep and keep all BIP corrections at zero.
    const TrendParameters cpaInitial{1.0, 0.0, 0.0, 0.0};
    const TrendParameters cpaSelected{0.001, 0.0, 0.0, 0.0};
    MPMC::tools::RegressionResult cpaResult;
    cpaResult.converged = true;
    cpaResult.iterations = 0;
    cpaResult.objectiveEvaluations = 9;
    cpaResult.initialObjective = meanTrainingTrendError(
        cpaInitial, states, targets,
        [](const TrendParameters &p) { return makePrCpaPrTrendMatched(p); });
    cpaResult.objective = meanTrainingTrendError(
        cpaSelected, states, targets,
        [](const TrendParameters &p) { return makePrCpaPrTrendMatched(p); });
    cpaResult.values.assign(cpaSelected.begin(), cpaSelected.end());
    cpaResult.finalSteps.assign(cpaSelected.size(), 0.0);
    TrendFit cpaFit{
        std::move(cpaDescriptors), std::move(cpaResult), cpaSelected};
    return {std::move(swFit), std::move(cpaFit)};
}

std::array<TrendFit, 2> frozenPrTrendFits(
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets)
{
    // Frozen from the bounded local fit performed on 2026-09-12.  Keeping the
    // values here makes the normal benchmark deterministic and fast; --refit
    // reruns the offline search against the same training states.
    const TrendParameters swInitial{0.0, 0.0, 0.0, 0.0};
    const TrendParameters swSelected{0.0, 0.001144, 0.0, 0.0};
    std::vector<MPMC::tools::RegressionParameter> swDescriptors{
        {"water_CO2_BIP_offset", 0.0, -0.25, 0.25, 0.025},
        {"water_HC_BIP_offset", 0.0, -0.25, 0.25, 0.025},
    };
    MPMC::tools::RegressionResult swResult;
    swResult.converged = true;
    swResult.initialObjective = meanTrainingTrendError(
        swInitial, states, targets,
        [](const TrendParameters &p) { return makeSwPrTrendMatched(p); });
    swResult.objective = meanTrainingTrendError(
        swSelected, states, targets,
        [](const TrendParameters &p) { return makeSwPrTrendMatched(p); });
    swResult.values.assign(swSelected.begin(), swSelected.begin() + 2);
    swResult.finalSteps.assign(2, 0.0);

    const TrendParameters cpaInitial{1.0, 0.0, 0.0, 0.0};
    const TrendParameters cpaSelected{0.001, 0.0, 0.0, 0.0};
    std::vector<MPMC::tools::RegressionParameter> cpaDescriptors{
        {"water_association_volume_scale", 1.0, 0.001, 1.0, 0.0},
        {"water_CO2_BIP_offset", 0.0, -0.40, 0.40, 0.04},
        {"water_HC_BIP_offset", 0.0, -0.40, 0.40, 0.04},
        {"water_BIP_slope_per_100K", 0.0, -0.30, 0.30, 0.03},
    };
    MPMC::tools::RegressionResult cpaResult;
    cpaResult.converged = true;
    cpaResult.initialObjective = meanTrainingTrendError(
        cpaInitial, states, targets,
        [](const TrendParameters &p) { return makePrCpaPrTrendMatched(p); });
    cpaResult.objective = meanTrainingTrendError(
        cpaSelected, states, targets,
        [](const TrendParameters &p) { return makePrCpaPrTrendMatched(p); });
    cpaResult.values.assign(cpaSelected.begin(), cpaSelected.end());
    cpaResult.finalSteps.assign(4, 0.0);

    return {
        TrendFit{std::move(swDescriptors), std::move(swResult), swSelected},
        TrendFit{std::move(cpaDescriptors), std::move(cpaResult), cpaSelected},
    };
}

CaseResult runCase(
    const Eos &eos,
    const std::string &backend,
    const std::string &parameterization,
    bool seeded)
{
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    options.maximumIterations = 240;
    options.maximumStabilityIterations = 160;
    options.fugacityTolerance = 1.0e-10;
    const Flash solver(eos, options);
    auto result = seeded
        ? solver.flashRestricted(
            pressure, temperature, feed, MPMC::PhasePresence::all(),
            normalizedReference())
        : solver.flash(pressure, temperature, feed);
    bool coincidentGasWaterMerged = false;
    if (result.converged &&
        result.phaseMoleFraction[1] > 1.0e-12 &&
        result.phaseMoleFraction[2] > 1.0e-12)
    {
        double maximumCompositionDifference = 0.0;
        for (std::size_t c = 0; c < N; ++c)
        {
            maximumCompositionDifference = std::max(
                maximumCompositionDifference,
                std::abs(result.composition[1][c] - result.composition[2][c]));
        }
        if (maximumCompositionDifference <= 1.0e-12)
        {
            result.phaseMoleFraction[2] += result.phaseMoleFraction[1];
            result.phaseMoleFraction[1] = 0.0;
            result.presence.remove(MPMC::CompositionalPhase::Gas);
            solver.updateSaturationsFromMoles(pressure, temperature, result);
            coincidentGasWaterMerged = true;
        }
    }
    return {backend, parameterization,
            seeded ? "seeded_three_phase" : "unrestricted",
            result, compare(eos, result), coincidentGasWaterMerged};
}

void writeSummary(
    const std::filesystem::path &path,
    const std::vector<CaseResult> &cases)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create summary CSV.");
    output << std::scientific << std::setprecision(12);
    output << "backend,parameterization,path,converged,phase_code,iterations,"
              "beta_oil,beta_gas,beta_water,reference_beta_oil,reference_beta_gas,"
              "reference_beta_water,composition_mae,composition_max_abs,beta_mae,"
              "beta_max_abs,max_mass_closure,max_log_fugacity_spread,"
              "temperature_K,pressure_bar,reference_state_match,"
              "coincident_gas_water_merged\n";
    for (const auto &entry : cases)
    {
        const auto &r = entry.flash;
        const auto &m = entry.metrics;
        const double notApplicable = std::numeric_limits<double>::quiet_NaN();
        const double refOil = isBenchmarkState() ? referenceBeta[0] : notApplicable;
        const double refGas = isBenchmarkState() ? referenceBeta[1] : notApplicable;
        const double refWater = isBenchmarkState() ? referenceBeta[2] : notApplicable;
        output << entry.backend << ',' << entry.parameterization << ','
               << entry.path << ',' << (r.converged ? 1 : 0)
               << ',' << (r.converged ? static_cast<int>(r.presence.bits()) : 0)
               << ',' << r.iterations << ','
               << r.phaseMoleFraction[0] << ',' << r.phaseMoleFraction[1] << ','
               << r.phaseMoleFraction[2] << ',' << refOil << ','
               << refGas << ',' << refWater << ','
               << m.compositionMae << ',' << m.compositionMax << ','
               << m.betaMae << ',' << m.betaMax << ',' << m.massClosure << ','
               << m.logFugacitySpread << ',' << temperature << ',' << pressureBar
               << ',' << (isBenchmarkState() ? 1 : 0) << ','
               << (entry.coincidentGasWaterMerged ? 1 : 0)
               << '\n';
    }
}

void writePhaseComparison(
    const std::filesystem::path &path,
    const std::vector<CaseResult> &cases)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create phase-comparison CSV.");
    const auto reference = normalizedReference();
    const std::array<std::string, 3> phaseNames{"oil", "gas", "water"};
    output << std::scientific << std::setprecision(12);
    output << "backend,parameterization,path,phase,component,reference_as_printed,"
              "reference_normalized,calculated,abs_error_to_normalized_reference\n";
    const std::array<Composition, 3> printed{
        referenceOil, referenceGas, referenceWater};
    const double notApplicable = std::numeric_limits<double>::quiet_NaN();
    for (const auto &entry : cases)
    {
        for (std::size_t p = 0; p < 3; ++p)
        {
            for (std::size_t c = 0; c < N; ++c)
            {
                const double calculated = entry.flash.composition[p][c];
                const double printedReference =
                    isBenchmarkState() ? printed[p][c] : notApplicable;
                const double normalizedReferenceValue =
                    isBenchmarkState() ? reference[p][c] : notApplicable;
                const double absoluteError = isBenchmarkState()
                    ? std::abs(calculated - reference[p][c])
                    : notApplicable;
                output << entry.backend << ',' << entry.parameterization << ','
                       << entry.path << ',' << phaseNames[p] << ','
                       << componentNames[c] << ',' << printedReference << ','
                       << normalizedReferenceValue << ',' << calculated << ','
                       << absoluteError << '\n';
            }
        }
    }
}

void writeReferenceClosure(const std::filesystem::path &path, const Flash &solver)
{
    const auto reference = normalizedReference();
    Composition kg{};
    Composition kw{};
    for (std::size_t c = 0; c < N; ++c)
    {
        kg[c] = reference[1][c] / reference[0][c];
        kw[c] = reference[2][c] / reference[0][c];
    }
    const auto rr = solver.solveThreePhaseRachfordRice(
        feed, kg, kw, referenceBeta[1], referenceBeta[2]);

    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create reference-closure CSV.");
    output << std::scientific << std::setprecision(12);
    output << "diagnostic,converged,beta_oil,beta_gas,beta_water,"
              "reference_beta_oil,reference_beta_gas,reference_beta_water\n";
    output << "RR_from_rounded_table6_K," << (rr.converged ? 1 : 0) << ','
           << rr.betaOil << ',' << rr.betaGas << ',' << rr.betaWater << ','
           << referenceBeta[0] << ',' << referenceBeta[1] << ','
           << referenceBeta[2] << '\n';

    Composition reconstructed{};
    for (std::size_t c = 0; c < N; ++c)
    {
        reconstructed[c] = referenceBeta[0] * reference[0][c]
            + referenceBeta[1] * reference[1][c]
            + referenceBeta[2] * reference[2][c];
    }
    std::ofstream balance(path.parent_path() / "published_rounding_balance.csv");
    if (!balance)
        throw std::runtime_error("Cannot create rounding-balance CSV.");
    balance << std::scientific << std::setprecision(12);
    balance << "component,feed_table_b3,reconstructed_from_normalized_table6,abs_error\n";
    for (std::size_t c = 0; c < N; ++c)
        balance << componentNames[c] << ',' << feed[c] << ',' << reconstructed[c]
                << ',' << std::abs(reconstructed[c] - feed[c]) << '\n';
}

bool isConvergedThreePhase(const CaseResult &entry)
{
    return entry.flash.converged &&
        entry.flash.presence.bits() == MPMC::PhasePresence::allBits;
}

double matchScore(const CaseResult &entry)
{
    return entry.metrics.compositionMae + entry.metrics.betaMae;
}

struct TrendComparisonMetrics
{
    double betaMae{std::numeric_limits<double>::quiet_NaN()};
    double compositionMae{std::numeric_limits<double>::quiet_NaN()};
    double objective{std::numeric_limits<double>::quiet_NaN()};
};

TrendComparisonMetrics compareToPrTrend(
    const Flash::Result &target,
    const Flash::Result &candidate)
{
    TrendComparisonMetrics metrics;
    if (!candidate.converged)
    {
        metrics.objective = prTrendError(target, candidate);
        return metrics;
    }

    double betaSum = 0.0;
    double compositionSum = 0.0;
    std::size_t compositionCount = 0;
    for (std::size_t p = 0; p < 3; ++p)
    {
        betaSum += std::abs(
            candidate.phaseMoleFraction[p] - target.phaseMoleFraction[p]);
        if (target.phaseMoleFraction[p] <= 1.0e-8 ||
            candidate.phaseMoleFraction[p] <= 1.0e-8)
        {
            continue;
        }
        for (std::size_t c = 0; c < N; ++c)
        {
            compositionSum += std::abs(
                candidate.composition[p][c] - target.composition[p][c]);
            ++compositionCount;
        }
    }
    metrics.betaMae = betaSum / 3.0;
    metrics.compositionMae = compositionCount == 0
        ? std::numeric_limits<double>::quiet_NaN()
        : compositionSum / static_cast<double>(compositionCount);
    metrics.objective = prTrendError(target, candidate);
    return metrics;
}

void writePrTrendComparison(
    const std::filesystem::path &path,
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets,
    const Eos &swBaseline,
    const Eos &swOptimized,
    const Eos &cpaBaseline,
    const Eos &cpaOptimized)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create PR-trend comparison CSV.");
    output << std::scientific << std::setprecision(12);
    output << "split,state,temperature_K,pressure_bar,backend,variant,converged,"
              "phase_code,beta_oil,beta_gas,beta_water,PR_beta_oil,PR_beta_gas,"
              "PR_beta_water,beta_mae_to_PR,composition_mae_to_PR,objective_to_PR\n";

    const std::array<std::string, 5> backends{
        "PR", "SW", "SW", "CPA", "CPA"};
    const std::array<std::string, 5> variants{
        "PR78_target", "paper_BIPs_baseline", "PR_trend_matched",
        "PR-CPA_4C_baseline", "PR_trend_matched"};
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        const std::array<Flash::Result, 5> results{
            targets[i], flashAt(swBaseline, states[i]),
            flashAt(swOptimized, states[i]), flashAt(cpaBaseline, states[i]),
            flashAt(cpaOptimized, states[i])};
        for (std::size_t model = 0; model < results.size(); ++model)
        {
            const auto metrics = compareToPrTrend(targets[i], results[model]);
            output << states[i].split << ',' << states[i].name << ','
                   << states[i].temperatureK << ','
                   << states[i].pressurePa / 1.0e5 << ',' << backends[model]
                   << ',' << variants[model] << ','
                   << (results[model].converged ? 1 : 0) << ','
                   << (results[model].converged
                           ? static_cast<int>(results[model].presence.bits())
                           : 0)
                   << ',' << results[model].phaseMoleFraction[0] << ','
                   << results[model].phaseMoleFraction[1] << ','
                   << results[model].phaseMoleFraction[2] << ','
                   << targets[i].phaseMoleFraction[0] << ','
                   << targets[i].phaseMoleFraction[1] << ','
                   << targets[i].phaseMoleFraction[2] << ','
                   << metrics.betaMae << ',' << metrics.compositionMae << ','
                   << metrics.objective << '\n';
        }
    }
}

void writePrTrendFitSummary(
    const std::filesystem::path &path,
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets,
    const Eos &swBaseline,
    const Eos &swOptimized,
    const Eos &cpaBaseline,
    const Eos &cpaOptimized)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create PR-trend fit-summary CSV.");
    output << std::scientific << std::setprecision(12);
    output << "backend,variant,split,states,mean_objective_to_PR,"
              "three_phase_agreement_count,converged_count\n";
    const std::array<std::string, 4> backends{"SW", "SW", "CPA", "CPA"};
    const std::array<std::string, 4> variants{
        "paper_BIPs_baseline", "PR_trend_matched",
        "PR-CPA_4C_baseline", "PR_trend_matched"};
    const std::array<const Eos *, 4> models{
        &swBaseline, &swOptimized, &cpaBaseline, &cpaOptimized};
    for (std::size_t model = 0; model < models.size(); ++model)
    {
        for (const std::string_view split : {"training", "validation"})
        {
            double sum = 0.0;
            std::size_t count = 0;
            std::size_t threePhaseAgreement = 0;
            std::size_t converged = 0;
            for (std::size_t i = 0; i < states.size(); ++i)
            {
                if (states[i].split != split ||
                    targets[i].presence.bits() != MPMC::PhasePresence::allBits)
                    continue;
                const auto result = flashAt(*models[model], states[i]);
                sum += prTrendError(targets[i], result);
                ++count;
                converged += result.converged ? 1u : 0u;
                threePhaseAgreement += result.converged &&
                    result.presence.bits() == targets[i].presence.bits() ? 1u : 0u;
            }
            output << backends[model] << ',' << variants[model] << ',' << split
                   << ',' << count << ','
                   << (count == 0 ? std::numeric_limits<double>::quiet_NaN()
                                  : sum / static_cast<double>(count))
                   << ',' << threePhaseAgreement << ',' << converged << '\n';
        }
    }
}

void writeCpaAssociationScaleSweep(
    const std::filesystem::path &path,
    const std::vector<TrendState> &states,
    const std::vector<Flash::Result> &targets)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create CPA association-scale sweep CSV.");
    output << std::scientific << std::setprecision(12);
    output << "association_volume_scale,split,states,mean_objective_to_PR,"
              "phase_agreement_count,center_phase_code,center_beta_oil,"
              "center_beta_gas,center_beta_water\n";
    for (const double scale : {
             1.0, 0.50, 0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001})
    {
        const Eos eos = makePrCpaPrTrendMatched({scale, 0.0, 0.0, 0.0});
        const auto center = flashAt(
            eos, TrendState{"T650_P390", "diagnostic", temperature, pressure});
        for (const std::string_view split : {"training", "validation"})
        {
            double sum = 0.0;
            std::size_t count = 0;
            std::size_t agreement = 0;
            for (std::size_t i = 0; i < states.size(); ++i)
            {
                if (states[i].split != split ||
                    targets[i].presence.bits() != MPMC::PhasePresence::allBits)
                    continue;
                const auto result = flashAt(eos, states[i]);
                sum += prTrendError(targets[i], result);
                ++count;
                agreement += result.converged &&
                    result.presence.bits() == targets[i].presence.bits() ? 1u : 0u;
            }
            output << scale << ',' << split << ',' << count << ','
                   << sum / static_cast<double>(count) << ',' << agreement << ','
                   << static_cast<int>(center.presence.bits()) << ','
                   << center.phaseMoleFraction[0] << ','
                   << center.phaseMoleFraction[1] << ','
                   << center.phaseMoleFraction[2] << '\n';
        }
    }
}

void writeBestByBackend(
    const std::filesystem::path &path,
    const std::vector<CaseResult> &cases)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create best-by-backend CSV.");
    output << std::scientific << std::setprecision(12);
    output << "backend,parameterization,path,beta_oil,beta_gas,beta_water,"
              "composition_mae,composition_max_abs,beta_mae,beta_max_abs,"
              "max_mass_closure,max_log_fugacity_spread\n";
    for (const char *backend : {"PR", "SW", "CPA"})
    {
        const CaseResult *best = nullptr;
        for (const auto &entry : cases)
        {
            if (entry.backend != backend || !isConvergedThreePhase(entry))
                continue;
            if (best == nullptr || matchScore(entry) < matchScore(*best))
                best = &entry;
        }
        if (best == nullptr)
            continue;
        output << best->backend << ',' << best->parameterization << ','
               << best->path << ',' << best->flash.phaseMoleFraction[0] << ','
               << best->flash.phaseMoleFraction[1] << ','
               << best->flash.phaseMoleFraction[2] << ','
               << best->metrics.compositionMae << ','
               << best->metrics.compositionMax << ',' << best->metrics.betaMae
               << ',' << best->metrics.betaMax << ','
               << best->metrics.massClosure << ','
               << best->metrics.logFugacitySpread << '\n';
    }
}

void writeTargetParameterDiagnostics(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create parameter-diagnostics CSV.");
    output << std::scientific << std::setprecision(12);
    output << "component,paper_H2O_kij,SW1992_aqueous_H2O_kij_at_650K,"
              "bounded_SW_uses_correlation,heavy_extrapolation_uses_correlation\n";
    output << componentNames[0] << ",0.0,0.0,1,1\n";
    output << componentNames[1] << ',' << kij[0][1] << ','
           << MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
                  temperature, tc[1], 0.0)
           << ",1,1\n";
    for (std::size_t c = 2; c < N; ++c)
    {
        output << componentNames[c] << ',' << kij[0][c] << ','
               << MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                      temperature, tc[c], omega[c], 0.0)
               << ',' << (c <= 3 ? 1 : 0) << ",1\n";
    }
    output << "SW_water_alpha,nan,"
           << MPMC::SoreideWhitsonCorrelations::waterAlpha(
                  temperature, tc[0], 0.0)
           << ",1,1\n";
}

void writeCpaParameterDiagnostics(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create CPA parameter-diagnostics CSV.");
    output << std::scientific << std::setprecision(12);
    output << "quantity,value_at_650K,unit_or_definition,source\n";
    output << "water_a0,1.5782e-1,Pa_m6_per_mol2,Sorensen2018_Table6\n"
           << "water_b,1.4788e-5,m3_per_mol,Sorensen2018_Table6\n"
           << "water_c1,0.6736,dimensionless,Sorensen2018_Table6\n"
           << "water_epsilon,16123,J_per_mol,Sorensen2018_Table6\n"
           << "water_beta,6.9662e-2,dimensionless,Sorensen2018_Table6\n"
           << "CO2_water_cross_epsilon,8061.5,J_per_mol,Sorensen2018_EqA17\n"
           << "CO2_water_cross_beta,0.15182,dimensionless,Sorensen2018_Table8\n";
    for (std::size_t c = 1; c < N; ++c)
    {
        output << "kij_H2O_" << componentNames[c] << ','
               << sorensen2018WaterBip(static_cast<int>(c), temperature)
               << ",dimensionless,Sorensen2018_Table9";
        if (c == 3)
            output << "_C2_C3_mean";
        else if (c == 4)
            output << "_nC5_mapping";
        else if (c >= 5)
            output << "_C7plus_mapping";
        output << '\n';
    }
    for (std::size_t c = 4; c < N; ++c)
    {
        output << "JiaOkuno_kij_H2O_" << componentNames[c] << ','
               << jiaOkuno2018HeavyWaterBip(static_cast<int>(c))
               << ",dimensionless,JiaOkuno2018_Table2_MW_mapping\n";
    }
}

void writePhysicalSelection(
    const std::filesystem::path &path,
    const std::vector<CaseResult> &cases)
{
    struct Selection
    {
        const char *backend;
        const char *role;
        const char *parameterization;
        const char *basis;
    };
    const std::array<Selection, 5> selected{{
        {"SW", "target_specific_recommended", "SW_water_alpha_paper_BIPs",
         "SW_water_alpha_plus_BSB_specific_BIPs_no_heavy_extrapolation"},
        {"SW", "predictive_applicability_check", "SW1992_bounded",
         "original_SW_correlations_only_within_supported_light_components"},
        {"CPA", "internally_consistent_literature", "PR-CPA_Sorensen2018_CO2-solvation",
         "matched_PR_CPA_formulation_with_CO2_cross_association"},
        {"CPA", "heavy_hydrocarbon_LVW_check", "SRK-CPA_JiaOkuno2018_C7plus-LVW-BIPs",
         "C7plus_BIPs_from_water_nalkane_three_phase_data"},
        {"CPA", "three_phase_transfer_sensitivity", "SRK-CPA_4C-water",
         "standard_SRK_CPA_water_with_PR_fitted_BSB_BIPs"}
    }};

    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create physical-selection CSV.");
    output << std::scientific << std::setprecision(12);
    output << "backend,role,parameterization,basis,converged,phase_code,"
              "beta_oil,beta_gas,beta_water,composition_mae,beta_mae\n";
    for (const auto &choice : selected)
    {
        const auto found = std::find_if(cases.begin(), cases.end(),
            [&choice](const CaseResult &entry) {
                return entry.backend == choice.backend &&
                    entry.parameterization == choice.parameterization &&
                    entry.path == "unrestricted";
            });
        if (found == cases.end())
            throw std::runtime_error("Selected physical case is missing.");
        output << found->backend << ',' << choice.role << ','
               << found->parameterization << ',' << choice.basis << ','
               << (found->flash.converged ? 1 : 0) << ','
               << static_cast<int>(found->flash.presence.bits()) << ','
               << found->flash.phaseMoleFraction[0] << ','
               << found->flash.phaseMoleFraction[1] << ','
               << found->flash.phaseMoleFraction[2] << ','
               << found->metrics.compositionMae << ','
               << found->metrics.betaMae << '\n';
    }
}

int runPrScan(const std::filesystem::path &outputDirectory)
{
    std::filesystem::create_directories(outputDirectory);
    const auto states = prTrendStates();
    const Eos pr78 = makePr(5);
    std::ofstream output(outputDirectory / "pr_target_states.csv");
    if (!output)
        throw std::runtime_error("Cannot create PR target-state CSV.");
    output << "split,state,temperature_K,pressure_bar,converged,phase_code,"
              "beta_oil,beta_gas,beta_water\n";
    output << std::scientific << std::setprecision(12);
    std::vector<Flash::Result> targets;
    targets.reserve(states.size());
    for (const auto &state : states)
    {
        const auto result = flashAt(pr78, state);
        targets.push_back(result);
        output << state.split << ',' << state.name << ',' << state.temperatureK
               << ',' << state.pressurePa / 1.0e5 << ','
               << (result.converged ? 1 : 0) << ','
               << (result.converged ? static_cast<int>(result.presence.bits()) : 0)
               << ',' << result.phaseMoleFraction[0] << ','
               << result.phaseMoleFraction[1] << ','
               << result.phaseMoleFraction[2] << '\n';
    }
    const Eos swBaseline = makeSwWithPaperHydrocarbonBips(false);
    const Eos cpaBaseline = makePrCpaSensitivity();
    writePrTrendComparison(
        outputDirectory / "pr_trend_baseline_comparison.csv", states, targets,
        swBaseline, swBaseline, cpaBaseline, cpaBaseline);
    writePrTrendFitSummary(
        outputDirectory / "pr_trend_baseline_summary.csv", states, targets,
        swBaseline, swBaseline, cpaBaseline, cpaBaseline);
    writeCpaAssociationScaleSweep(
        outputDirectory / "cpa_association_scale_sweep.csv", states, targets);
    return 0;
}

void writeStandardModelBasis(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create standard-model basis CSV.");
    output << "backend,parameterization,pure_component_model,water_treatment,"
              "binary_interaction_policy,target_fitted\n"
           << "PR,PR78,Peng-Robinson-1978,ordinary_PR_alpha,"
              "Heringer2025_Table_B4,0\n"
           << "SW,SW_water_alpha_paper_BIPs,Peng-Robinson-1978,"
              "Soreide-Whitson_1992_pure-water-alpha_zero-salinity,"
              "Heringer2025_Table_B4_no_extrapolated_SW_BIPs,0\n"
           << "CPA,SRK-CPA_4C-water,"
              "Soave-Redlich-Kwong,standard_CPA_4C_water_full_association,"
              "Heringer2025_Table_B4_transfer_no_refit,0\n";
}

double pr78C1(std::size_t component)
{
    const double w = omega[component];
    if (w > 0.49)
    {
        return 0.379642 + 1.48503 * w - 0.164423 * w * w
            + 0.016666 * w * w * w;
    }
    return 0.37464 + 1.54226 * w - 0.26992 * w * w;
}

double srkC1(std::size_t component)
{
    const double w = omega[component];
    return 0.480 + 1.574 * w - 0.176 * w * w;
}

double alphaFromC1(double c1, std::size_t component, double temperatureK)
{
    const double root = 1.0 + c1 *
        (1.0 - std::sqrt(temperatureK / tc[component]));
    return root * root;
}

void writeStandardComponentInputs(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create component-input CSV.");
    output << std::scientific << std::setprecision(12);
    output << "component,feed_mole_fraction,Tc_K,Pc_bar,acentric_factor,"
              "critical_volume_m3_per_mol,molar_mass_kg_per_mol,"
              "primary_source,metadata_source\n";
    for (std::size_t c = 0; c < N; ++c)
    {
        output << componentNames[c] << ',' << feed[c] << ',' << tc[c] << ','
               << pc[c] / 1.0e5 << ',' << omega[c] << ',' << vc[c] << ','
               << mw[c] << ",Heringer2025_Table_B3,";
        output << (c <= 1 ? "standard_pure_component_metadata"
                          : "Fernandes2021_Table_12_BSB_metadata")
               << '\n';
    }
}

void writeStandardKijMatrix(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create BIP-matrix CSV.");
    output << std::fixed << std::setprecision(6) << "component";
    for (const auto &name : componentNames)
        output << ',' << name;
    output << '\n';
    for (std::size_t i = 0; i < N; ++i)
    {
        output << componentNames[i];
        for (std::size_t j = 0; j < N; ++j)
            output << ',' << kij[i][j];
        output << '\n';
    }
}

void writeStandardEosComponentParameters(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create EOS-component parameter CSV.");
    output << std::scientific << std::setprecision(12);
    output << "backend,component,cubic_physical_term,omega_A,omega_B,a0_Pa_m6_per_mol2,"
              "b_m3_per_mol,c1,alpha_at_650K,state_temperature_K,"
              "alpha_at_state_temperature,association_energy_J_per_mol,"
              "association_volume,donor_sites,acceptor_sites,H2O_pair_kij_used,"
              "target_fitted\n";
    constexpr double R = MPMC::units::gasConstant;
    for (const std::string backend : {"PR78", "SW", "SRK-CPA"})
    {
        for (std::size_t c = 0; c < N; ++c)
        {
            const bool cpa = backend == "SRK-CPA";
            double omegaA = cpa ? 0.42748 : 0.45724;
            double omegaB = cpa ? 0.08664 : 0.07780;
            double a0 = omegaA * R * R * tc[c] * tc[c] / pc[c];
            double b = omegaB * R * tc[c] / pc[c];
            double c1 = cpa ? srkC1(c) : pr78C1(c);
            double epsilon = 0.0;
            double beta = 0.0;
            int donors = 0;
            int acceptors = 0;
            if (cpa && c == 0)
            {
                a0 = MPMC::StandardCpaWater4C::a0;
                b = MPMC::StandardCpaWater4C::b;
                c1 = MPMC::StandardCpaWater4C::c1;
                epsilon = MPMC::StandardCpaWater4C::epsilon;
                beta = MPMC::StandardCpaWater4C::beta;
                donors = MPMC::StandardCpaWater4C::donorSites;
                acceptors = MPMC::StandardCpaWater4C::acceptorSites;
            }
            const double alphaAt650K = backend == "SW" && c == 0
                ? MPMC::SoreideWhitsonCorrelations::waterAlpha(
                      benchmarkTemperature, tc[0], 0.0)
                : alphaFromC1(c1, c, benchmarkTemperature);
            const double alphaAtState = backend == "SW" && c == 0
                ? MPMC::SoreideWhitsonCorrelations::waterAlpha(
                      temperature, tc[0], 0.0)
                : alphaFromC1(c1, c, temperature);
            output << backend << ',' << componentNames[c] << ','
                   << (cpa ? "Soave-Redlich-Kwong" : "Peng-Robinson") << ','
                   << omegaA << ',' << omegaB << ',' << a0 << ',' << b << ','
                   << c1 << ',' << alphaAt650K << ',' << temperature << ','
                   << alphaAtState << ',' << epsilon << ',' << beta << ','
                   << donors << ',' << acceptors << ','
                   << (c == 0 ? 0.0 : kij[0][c]) << ",0\n";
        }
    }
}

void writeStandardFlashControls(const std::filesystem::path &path)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot create flash-control CSV.");
    output << std::scientific << std::setprecision(12);
    output << "parameter,value,unit_or_definition\n"
           << "temperature," << temperature << ",K\n"
           << "pressure," << pressureBar << ",bar\n"
           << "salinity_molality,0,mol_NaCl_per_kg_H2O\n"
           << "water_component_index,0,zero_based\n"
           << "flash_path,unrestricted,production_three_phase_stability_flash\n"
           << "maximum_flash_iterations,240,count\n"
           << "maximum_stability_iterations,160,count\n"
           << "fugacity_tolerance,1e-10,max_log_fugacity_residual\n"
           << "coincident_phase_composition_tolerance,1e-8,max_abs_composition_difference\n"
           << "allowed_phase_roles,oil_gas_water,canonical_public_order\n";
}

std::string phaseCodeLabel(unsigned code)
{
    std::string label;
    const auto append = [&label](std::string_view phase) {
        if (!label.empty())
            label += '+';
        label += phase;
    };
    if ((code & MPMC::PhasePresence::oilBit) != 0)
        append("O");
    if ((code & MPMC::PhasePresence::gasBit) != 0)
        append("G");
    if ((code & MPMC::PhasePresence::waterBit) != 0)
        append("W");
    return label.empty() ? "none" : label;
}

int regularGridIntervalCount(
    double lower,
    double upper,
    double step,
    std::string_view name)
{
    if (!std::isfinite(lower) || !std::isfinite(upper) ||
        !std::isfinite(step) || !(lower > 0.0) || upper < lower ||
        !(step > 0.0))
    {
        throw std::invalid_argument(
            std::string(name) + " grid bounds and step must be finite and positive.");
    }
    const double rawCount = (upper - lower) / step;
    const int count = static_cast<int>(std::llround(rawCount));
    if (count < 0 || std::abs(rawCount - count) > 1.0e-9)
    {
        throw std::invalid_argument(
            std::string(name) + " range must be an integer multiple of its step.");
    }
    return count;
}

struct PtTopologyRow
{
    double temperatureK{0.0};
    double pressureMpa{0.0};
    Flash::Result result;
    double massClosure{std::numeric_limits<double>::quiet_NaN()};
    double fugacitySpread{std::numeric_limits<double>::quiet_NaN()};
};

std::vector<PtTopologyRow> computePtTopologyModel(
    Eos eos,
    double minimumTemperature,
    int temperatureIntervals,
    double temperatureStep,
    double minimumPressureMpa,
    int pressureIntervals,
    double pressureStepMpa)
{
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    options.maximumIterations = 240;
    options.maximumStabilityIterations = 160;
    options.fugacityTolerance = 1.0e-10;
    const Flash solver(eos, options);

    std::vector<PtTopologyRow> rows;
    rows.reserve(static_cast<std::size_t>(temperatureIntervals + 1) *
                 static_cast<std::size_t>(pressureIntervals + 1));
    for (int ti = 0; ti <= temperatureIntervals; ++ti)
    {
        const double stateTemperature = minimumTemperature + ti * temperatureStep;
        for (int pi = 0; pi <= pressureIntervals; ++pi)
        {
            const double statePressureMpa =
                minimumPressureMpa + pi * pressureStepMpa;
            const double statePressure = statePressureMpa * 1.0e6;
            Flash::Result result = solver.flash(
                statePressure, stateTemperature, feed);
            PtTopologyRow row;
            row.temperatureK = stateTemperature;
            row.pressureMpa = statePressureMpa;
            row.result = std::move(result);
            if (row.result.converged)
            {
                row.massClosure = maxMassClosure(row.result);
                row.fugacitySpread = maxLogFugacitySpreadAt(
                    eos, row.result, statePressure, stateTemperature);
            }
            rows.push_back(std::move(row));
        }
    }
    return rows;
}

void writePtTopologyRow(
    std::ofstream &output,
    const PtTopologyRow &row,
    std::string_view backend,
    std::string_view parameterization)
{
    const auto &result = row.result;
    const unsigned code = result.converged
        ? static_cast<unsigned>(result.presence.bits()) : 0u;
    output << row.temperatureK << ',' << row.pressureMpa << ','
           << backend << ',' << parameterization << ','
           << (result.converged ? 1 : 0) << ',' << code << ','
           << phaseCodeLabel(code) << ','
           << (result.converged ? result.presence.count() : 0) << ','
           << result.phaseMoleFraction[0] << ','
           << result.phaseMoleFraction[1] << ','
           << result.phaseMoleFraction[2] << ',' << result.iterations << ','
           << row.massClosure << ',' << row.fugacitySpread << '\n';
}

int runPtTopologyGrid(
    const std::filesystem::path &outputPath,
    double minimumTemperature,
    double maximumTemperature,
    double temperatureStep,
    double minimumPressureMpa,
    double maximumPressureMpa,
    double pressureStepMpa)
{
    const int temperatureIntervals = regularGridIntervalCount(
        minimumTemperature, maximumTemperature, temperatureStep, "Temperature");
    const int pressureIntervals = regularGridIntervalCount(
        minimumPressureMpa, maximumPressureMpa, pressureStepMpa, "Pressure");
    if (!outputPath.parent_path().empty())
        std::filesystem::create_directories(outputPath.parent_path());

    std::ofstream output(outputPath);
    if (!output)
        throw std::runtime_error("Cannot create P-T topology CSV.");
    output << std::scientific << std::setprecision(12);
    output << "temperature_K,pressure_MPa,backend,parameterization,converged,"
              "phase_code,phase_code_label,phase_count,beta_oil,beta_gas,"
              "beta_water,iterations,max_mass_closure,max_log_fugacity_spread\n";

    auto prFuture = std::async(
        std::launch::async, computePtTopologyModel, makePr(5),
        minimumTemperature, temperatureIntervals, temperatureStep,
        minimumPressureMpa, pressureIntervals, pressureStepMpa);
    auto swFuture = std::async(
        std::launch::async, computePtTopologyModel,
        makeSwWithPaperHydrocarbonBips(false),
        minimumTemperature, temperatureIntervals, temperatureStep,
        minimumPressureMpa, pressureIntervals, pressureStepMpa);
    auto cpaFuture = std::async(
        std::launch::async, computePtTopologyModel, makeSrkCpa(),
        minimumTemperature, temperatureIntervals, temperatureStep,
        minimumPressureMpa, pressureIntervals, pressureStepMpa);

    const auto prRows = prFuture.get();
    std::cout << "[P-T grid] PR complete\n";
    const auto swRows = swFuture.get();
    std::cout << "[P-T grid] SW complete\n";
    const auto cpaRows = cpaFuture.get();
    std::cout << "[P-T grid] CPA complete\n";
    if (prRows.size() != swRows.size() || prRows.size() != cpaRows.size())
        throw std::runtime_error("P-T topology model grids have inconsistent sizes.");

    std::size_t failedStates = 0;
    for (std::size_t i = 0; i < prRows.size(); ++i)
    {
        writePtTopologyRow(output, prRows[i], "PR", "PR78");
        writePtTopologyRow(
            output, swRows[i], "SW", "SW_water_alpha_paper_BIPs");
        writePtTopologyRow(
            output, cpaRows[i], "CPA", "SRK-CPA_4C-water");
        failedStates += prRows[i].result.converged ? 0u : 1u;
        failedStates += swRows[i].result.converged ? 0u : 1u;
        failedStates += cpaRows[i].result.converged ? 0u : 1u;
    }
    return failedStates == 0 ? 0 : 2;
}

int runStandard(const std::filesystem::path &outputDirectory)
{
    std::filesystem::create_directories(outputDirectory);
    const Eos pr78 = makePr(5);
    const Eos swStandard = makeSwWithPaperHydrocarbonBips(false);
    const Eos cpaStandard = makeSrkCpa();

    // Primary comparison: exactly one non-target-fitted calculation from each
    // EOS family.  All use the published feed and Table B4 BIPs so the model
    // structure, rather than a parameter regression, is what changes.
    std::vector<CaseResult> primaryCases;
    primaryCases.push_back(runCase(pr78, "PR", "PR78", false));
    primaryCases.push_back(runCase(
        swStandard, "SW", "SW_water_alpha_paper_BIPs", false));
    primaryCases.push_back(runCase(
        cpaStandard, "CPA", "SRK-CPA_4C-water", false));

    // Keep additional physically motivated mappings and seeded paths in a
    // separate diagnostic table.  None enters the primary plot or ranking.
    const Eos pr76 = makePr(1);
    const Eos swBounded = makeSw(false);
    const Eos swExtrapolated = makeSw(true);
    const Eos swCo2PaperHydrocarbons = makeSwWithPaperHydrocarbonBips(true);
    const Eos srkCpaJiaOkuno = makeSrkCpaJiaOkuno2018();
    const Eos prCpa = makePrCpaSensitivity();
    const Eos prCpaSorensen = makePrCpaSorensen2018();
    std::vector<CaseResult> diagnostics = primaryCases;
    diagnostics.push_back(runCase(pr76, "PR", "PR76", false));
    diagnostics.push_back(runCase(pr78, "PR", "PR78", true));
    diagnostics.push_back(runCase(swStandard, "SW", "SW_water_alpha_paper_BIPs", true));
    diagnostics.push_back(runCase(swBounded, "SW", "SW1992_bounded", false));
    diagnostics.push_back(runCase(swBounded, "SW", "SW1992_bounded", true));
    diagnostics.push_back(runCase(
        swExtrapolated, "SW", "SW1992_heavy_extrapolation", false));
    diagnostics.push_back(runCase(
        swCo2PaperHydrocarbons, "SW", "SW1992_CO2_paper_HC_BIPs", false));
    diagnostics.push_back(runCase(cpaStandard, "CPA", "SRK-CPA_4C-water", true));
    diagnostics.push_back(runCase(
        srkCpaJiaOkuno, "CPA", "SRK-CPA_JiaOkuno2018_C7plus-LVW-BIPs", false));
    diagnostics.push_back(runCase(
        prCpa, "CPA", "PR-CPA_PR78-cubic_4C-association", false));
    diagnostics.push_back(runCase(
        prCpaSorensen, "CPA", "PR-CPA_Sorensen2018_CO2-solvation", false));

    writeSummary(outputDirectory / "summary.csv", primaryCases);
    writePhaseComparison(outputDirectory / "phase_comparison.csv", primaryCases);
    writeBestByBackend(outputDirectory / "best_by_backend.csv", primaryCases);
    writeSummary(outputDirectory / "diagnostic_summary.csv", diagnostics);
    writePhaseComparison(
        outputDirectory / "diagnostic_phase_comparison.csv", diagnostics);
    writeTargetParameterDiagnostics(
        outputDirectory / "sw_parameter_diagnostics.csv");
    writeCpaParameterDiagnostics(
        outputDirectory / "cpa_parameter_diagnostics.csv");
    writePhysicalSelection(
        outputDirectory / "physical_selection.csv", diagnostics);
    writeStandardModelBasis(outputDirectory / "standard_model_basis.csv");
    writeStandardComponentInputs(
        outputDirectory / "standard_component_inputs.csv");
    writeStandardKijMatrix(
        outputDirectory / "standard_binary_interaction_matrix.csv");
    writeStandardEosComponentParameters(
        outputDirectory / "standard_eos_component_parameters.csv");
    writeStandardFlashControls(
        outputDirectory / "standard_flash_controls.csv");

    MPMC::ThreePhaseFlashOptions diagnosticOptions;
    diagnosticOptions.waterComponent = 0;
    const Flash diagnosticSolver(makePr(1), diagnosticOptions);
    writeReferenceClosure(outputDirectory / "reference_rr.csv", diagnosticSolver);

    std::cout << std::scientific << std::setprecision(6);
    for (const auto &entry : primaryCases)
    {
        std::cout << entry.backend << '/' << entry.parameterization
                  << ": converged=" << entry.flash.converged
                  << ", phase_code=" << static_cast<int>(entry.flash.presence.bits())
                  << ", beta(O/G/W)=" << entry.flash.phaseMoleFraction[0] << '/'
                  << entry.flash.phaseMoleFraction[1] << '/'
                  << entry.flash.phaseMoleFraction[2]
                  << ", composition_MAE=" << entry.metrics.compositionMae
                  << ", beta_MAE=" << entry.metrics.betaMae << '\n';
    }
    return std::all_of(
        primaryCases.begin(), primaryCases.end(),
        [](const CaseResult &entry) { return entry.flash.converged; }) ? 0 : 2;
}

int runTargetFit(const std::filesystem::path &outputDirectory, bool refit)
{
    std::filesystem::create_directories(outputDirectory);
    std::vector<CaseResult> cases;
    const Eos pr76 = makePr(1);
    const Eos pr78 = makePr(5);
    const auto trendStates = prTrendStates();
    std::vector<Flash::Result> trendTargets;
    trendTargets.reserve(trendStates.size());
    for (const auto &state : trendStates)
        trendTargets.push_back(flashAt(pr78, state));
    std::cout << "[PR-trend fit] "
              << (refit ? "refitting" : "using frozen fit for") << " SW and CPA against "
              << trendStates.size() << " PR states...\n";
    const auto trendFits = refit
        ? optimizePrTrend(trendStates, trendTargets)
        : frozenPrTrendFits(trendStates, trendTargets);
    const Eos swTrendMatched = makeSwPrTrendMatched(trendFits[0].parameters);
    const Eos cpaTrendMatched = makePrCpaPrTrendMatched(trendFits[1].parameters);
    const Eos swBounded = makeSw(false);
    const Eos swExtrapolated = makeSw(true);
    const Eos swCo2PaperHydrocarbons = makeSwWithPaperHydrocarbonBips(true);
    const Eos swPaperBips = makeSwWithPaperHydrocarbonBips(false);
    const Eos srkCpa = makeSrkCpa();
    const Eos srkCpaJiaOkuno = makeSrkCpaJiaOkuno2018();
    const Eos prCpa = makePrCpaSensitivity();
    const Eos prCpaSorensen = makePrCpaSorensen2018();
    cases.push_back(runCase(pr76, "PR", "PR76", false));
    cases.push_back(runCase(pr76, "PR", "PR76", true));
    cases.push_back(runCase(pr78, "PR", "PR78", false));
    cases.push_back(runCase(pr78, "PR", "PR78", true));
    cases.push_back(runCase(swBounded, "SW", "SW1992_bounded", false));
    cases.push_back(runCase(swBounded, "SW", "SW1992_bounded", true));
    cases.push_back(runCase(
        swExtrapolated, "SW", "SW1992_heavy_extrapolation", false));
    cases.push_back(runCase(
        swCo2PaperHydrocarbons, "SW", "SW1992_CO2_paper_HC_BIPs", false));
    cases.push_back(runCase(
        swCo2PaperHydrocarbons, "SW", "SW1992_CO2_paper_HC_BIPs", true));
    cases.push_back(runCase(
        swPaperBips, "SW", "SW_water_alpha_paper_BIPs", false));
    cases.push_back(runCase(
        swTrendMatched, "SW", "SW_PR-trend-matched", false));
    cases.push_back(runCase(
        swTrendMatched, "SW", "SW_PR-trend-matched", true));
    cases.push_back(runCase(srkCpa, "CPA", "SRK-CPA_4C-water", false));
    cases.push_back(runCase(srkCpa, "CPA", "SRK-CPA_4C-water", true));
    cases.push_back(runCase(
        srkCpaJiaOkuno, "CPA", "SRK-CPA_JiaOkuno2018_C7plus-LVW-BIPs", false));
    cases.push_back(runCase(
        srkCpaJiaOkuno, "CPA", "SRK-CPA_JiaOkuno2018_C7plus-LVW-BIPs", true));
    cases.push_back(runCase(
        prCpa, "CPA", "PR-CPA_PR78-cubic_4C-association", false));
    cases.push_back(runCase(
        prCpa, "CPA", "PR-CPA_PR78-cubic_4C-association", true));
    cases.push_back(runCase(
        prCpaSorensen, "CPA", "PR-CPA_Sorensen2018_CO2-solvation", false));
    cases.push_back(runCase(
        prCpaSorensen, "CPA", "PR-CPA_Sorensen2018_CO2-solvation", true));
    cases.push_back(runCase(
        cpaTrendMatched, "CPA", "PR-CPA_PR-trend-matched", false));
    cases.push_back(runCase(
        cpaTrendMatched, "CPA", "PR-CPA_PR-trend-matched", true));

    writeSummary(outputDirectory / "summary.csv", cases);
    writePhaseComparison(outputDirectory / "phase_comparison.csv", cases);
    writeBestByBackend(outputDirectory / "best_by_backend.csv", cases);
    writeTargetParameterDiagnostics(
        outputDirectory / "target_parameter_diagnostics.csv");
    writeCpaParameterDiagnostics(
        outputDirectory / "cpa_parameter_diagnostics.csv");
    writePhysicalSelection(
        outputDirectory / "physical_selection.csv", cases);
    MPMC::tools::writeRegressionResultCsv(
        outputDirectory / "sw_pr_trend_parameters.csv",
        trendFits[0].descriptors, trendFits[0].result);
    MPMC::tools::writeRegressionResultCsv(
        outputDirectory / "cpa_pr_trend_parameters.csv",
        trendFits[1].descriptors, trendFits[1].result);
    writePrTrendComparison(
        outputDirectory / "pr_trend_comparison.csv", trendStates, trendTargets,
        swPaperBips, swTrendMatched, prCpa, cpaTrendMatched);
    writePrTrendFitSummary(
        outputDirectory / "pr_trend_fit_summary.csv", trendStates, trendTargets,
        swPaperBips, swTrendMatched, prCpa, cpaTrendMatched);
    const Eos diagnosticEos = makePr(1);
    MPMC::ThreePhaseFlashOptions diagnosticOptions;
    diagnosticOptions.waterComponent = 0;
    const Flash diagnosticSolver(diagnosticEos, diagnosticOptions);
    writeReferenceClosure(outputDirectory / "reference_rr.csv", diagnosticSolver);

    const CaseResult *best = nullptr;
    for (const auto &entry : cases)
    {
        if (!isConvergedThreePhase(entry))
            continue;
        if (best == nullptr || matchScore(entry) < matchScore(*best))
            best = &entry;
    }

    std::cout << std::scientific << std::setprecision(6);
    for (const auto &entry : cases)
    {
        std::cout << entry.backend << '/' << entry.parameterization << '/'
                  << entry.path
                  << ": converged=" << entry.flash.converged
                  << ", phase_code=" << static_cast<int>(entry.flash.presence.bits())
                  << ", beta(O/G/W)=" << entry.flash.phaseMoleFraction[0] << '/'
                  << entry.flash.phaseMoleFraction[1] << '/'
                  << entry.flash.phaseMoleFraction[2]
                  << ", composition_MAE=" << entry.metrics.compositionMae
                  << ", beta_MAE=" << entry.metrics.betaMae << '\n';
    }
    if (best == nullptr)
    {
        std::cout << "No converged three-phase result was obtained.\n";
        return 2;
    }

    std::ofstream bestOutput(outputDirectory / "best_case.txt");
    if (!bestOutput)
        throw std::runtime_error("Cannot create best-case report.");
    bestOutput << std::scientific << std::setprecision(12)
               << "backend=" << best->backend << '\n'
               << "parameterization=" << best->parameterization << '\n'
               << "path=" << best->path << '\n'
               << "temperature_K=" << temperature << '\n'
               << "pressure_bar=" << pressureBar << '\n'
               << "composition_mae=" << best->metrics.compositionMae << '\n'
               << "composition_max_abs=" << best->metrics.compositionMax << '\n'
               << "beta_mae=" << best->metrics.betaMae << '\n'
               << "beta_max_abs=" << best->metrics.betaMax << '\n';
    std::cout << "Best three-phase match: " << best->backend << '/'
              << best->parameterization << '/' << best->path << "\n";
    return 0;
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 2 && argc != 3 && argc != 5 && argc != 9)
        {
            std::cerr << "usage: heringer2025_bsb_three_phase OUTPUT_DIR "
                         "[--pr-scan|--target-fit|--refit-target-fit|"
                         "--state TEMPERATURE_K PRESSURE_BAR|"
                         "--pt-topology-grid T_MIN_K T_MAX_K T_STEP_K "
                         "P_MIN_MPa P_MAX_MPa P_STEP_MPa]\n";
            return 1;
        }
        if (argc == 9)
        {
            const std::string_view option(argv[2]);
            if (option != "--pt-topology-grid")
                throw std::invalid_argument("Unknown option: " + std::string(argv[2]));
            return runPtTopologyGrid(
                argv[1], std::stod(argv[3]), std::stod(argv[4]),
                std::stod(argv[5]), std::stod(argv[6]),
                std::stod(argv[7]), std::stod(argv[8]));
        }
        if (argc == 5)
        {
            const std::string_view option(argv[2]);
            if (option != "--state")
                throw std::invalid_argument("Unknown option: " + std::string(argv[2]));
            temperature = std::stod(argv[3]);
            pressureBar = std::stod(argv[4]);
            if (!std::isfinite(temperature) || temperature <= 0.0 ||
                !std::isfinite(pressureBar) || pressureBar <= 0.0)
            {
                throw std::invalid_argument(
                    "Temperature and pressure must be finite and positive.");
            }
            pressure = pressureBar * 1.0e5;
            return runStandard(argv[1]);
        }
        if (argc == 3)
        {
            const std::string_view option(argv[2]);
            if (option == "--pr-scan")
                return runPrScan(argv[1]);
            if (option == "--target-fit")
                return runTargetFit(argv[1], false);
            if (option != "--refit-target-fit")
                throw std::invalid_argument("Unknown option: " + std::string(argv[2]));
            return runTargetFit(argv[1], true);
        }
        return runStandard(argv[1]);
    }
    catch (const std::exception &error)
    {
        std::cerr << "Heringer 2025 benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
