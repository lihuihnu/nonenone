/**
 * @file thermo_empirical_parameters_test.cpp
 * @brief 单元测试：验证 `thermo_empirical_parameters` 的核心语义、边界条件和回归行为。
 */
#include <common/units.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/aqueous_volume.hpp>
#include <natural/properties/aqueous_viscosity.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

using ModelConfig = MPMC::CompositionalModelConfig<
    2, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<ModelConfig>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Composition = std::array<double, 2>;
using AdIndices = MPMC::ADIndices<ModelConfig>;
using AdEos = MPMC::CubicEquationOfState<AdIndices>;
using AdValue = AdIndices::ValueType;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(message);
}

struct ExtendedFactoryConfig final
{
    struct Fluid
    {
        static constexpr std::size_t N = 2;
        inline static constexpr std::array<const char *, N> componentNames{"C2", "nC5"};
        inline static constexpr std::array<double, N> criticalTemperature{305.32, 469.7};
        inline static constexpr std::array<double, N> criticalPressure{4.872e6, 3.37e6};
        inline static constexpr std::array<double, N> criticalVolume{1.45e-4, 3.13e-4};
        inline static constexpr std::array<double, N> acentricFactor{0.0995, 0.251};
        inline static constexpr std::array<double, N> molarMass{0.03007, 0.07215};
        inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
            {{0.0, 0.0200}}, {{0.0200, 0.0}}
        }};
        inline static constexpr std::array<std::array<double, N>, N> binaryInteractionTemperatureSlope{{
            {{0.0, 2.0e-4}}, {{2.0e-4, 0.0}}
        }};
        static constexpr double binaryInteractionReferenceTemperature = 300.0;
        inline static constexpr std::array<double, N> componentVolumeTranslation{
            1.0e-6, 3.0e-6};

        static constexpr double eosOmegaA = 0.45724;
        static constexpr double eosOmegaB = 0.07780;
        static constexpr int eosModelFlag = 1;
        static constexpr double eosU = 2.414213562373095;
        static constexpr double eosW = -0.414213562373095;
        static constexpr double temperature = 340.0;
        static constexpr int waterComponent = 0;
        inline static constexpr std::array<double, 3> surfaceDensity{600.0, 30.0, 1000.0};
        inline static constexpr std::array<double, 3> viscosity{2.0e-4, 1.5e-5, 5.0e-4};
    };
};

struct CallbackFactoryConfig final
{
    struct Fluid : ExtendedFactoryConfig::Fluid
    {
        static double cubicBinaryInteractionCoefficient(
            int i, int j, double temperature)
        {
            if (i == j)
                return 0.0;
            return 0.125 + 2.5e-4 * (temperature - 310.0);
        }
    };
};

template <class TestIndices>
MPMC::CubicEquationOfState<TestIndices> makeB2Cpa()
{
    constexpr std::array<double, 2> tc{647.096, 304.1282};
    constexpr std::array<double, 2> pc{220.64e5, 73.773e5};
    constexpr std::array<double, 2> vc{5.6e-5, 9.4e-5};
    constexpr std::array<double, 2> omega{0.344, 0.22394};
    constexpr std::array<double, 2> mw{0.01801528, 0.0440095};
    constexpr std::array<std::array<double, 2>, 2> kij{{
        {{0.0, 0.11406}}, {{0.11406, 0.0}}
    }};
    MPMC::CubicEquationOfState<TestIndices> eos(
        0.42748, 0.08664,
        MPMC::CompositionalMixture<TestIndices>(tc, pc, vc, omega, mw, kij),
        5, 1.0, 0.0, 1.0e-30);
    typename MPMC::CubicEquationOfState<TestIndices>::CubicPlusAssociationOptions cpa;
    constexpr double R = MPMC::units::gasConstant;
    cpa.a0 = {14.52e-6 * R * 1017.3, 27.2e-6 * R * 1551.22};
    cpa.b = {14.52e-6, 27.2e-6};
    cpa.c1 = {0.6736, 0.7602};
    cpa.associationEnergy = {2003.2 * R, 14200.0};
    cpa.associationVolume = {0.0692, 0.0162};
    cpa.donorSites = {2, 0};
    cpa.acceptorSites = {2, 1};
    cpa.crossAssociationEnergy[0][1] = 14200.0;
    cpa.crossAssociationVolume[0][1] = 0.0162;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    eos.configureCpaTemperatureCache(305.0);
    return eos;
}

template <class TestIndices>
MPMC::CubicEquationOfState<TestIndices> makeWaterMethanolCpa()
{
    constexpr std::array<double, 2> tc{647.3, 512.6};
    constexpr std::array<double, 2> pc{22.0483e6, 8.0959e6};
    constexpr std::array<double, 2> vc{5.60e-5, 1.18e-4};
    constexpr std::array<double, 2> omega{0.344, 0.559};
    constexpr std::array<double, 2> mw{0.018015, 0.032042};
    constexpr std::array<std::array<double, 2>, 2> kij{{
        {{0.0, -0.09}}, {{-0.09, 0.0}}
    }};
    MPMC::CubicEquationOfState<TestIndices> eos(
        0.42748, 0.08664,
        MPMC::CompositionalMixture<TestIndices>(tc, pc, vc, omega, mw, kij),
        1, 1.0, 0.0, 1.0e-30);
    typename MPMC::CubicEquationOfState<TestIndices>::CubicPlusAssociationOptions cpa;
    cpa.a0 = {0.12277, 0.40531};
    cpa.b = {1.4515e-5, 3.0978e-5};
    cpa.c1 = {0.67359, 0.43102};
    cpa.associationEnergy = {16655.0, 24591.0};
    cpa.associationVolume = {0.0692, 0.0161};
    cpa.donorSites = {2, 1};
    cpa.acceptorSites = {2, 1};
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(std::move(cpa));
    eos.configureCpaTemperatureCache(389.375);
    return eos;
}

void checkFactoryTemperatureDependentBip()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, ExtendedFactoryConfig>();
    require(fluid.eos.usesTemperatureDependentBinaryInteraction(),
            "factory must enable the optional temperature-dependent BIP");

    near(fluid.eos.binaryInteractionCoefficient(0, 1, 300.0), 0.0200, 1.0e-15,
         "temperature-dependent BIP reference value");
    near(fluid.eos.binaryInteractionCoefficient(0, 1, 350.0), 0.0300, 1.0e-15,
         "temperature-dependent BIP linear slope");

    const auto at300 = fluid.eos.mixingParameters(
        5.0e6, 300.0, MPMC::CompositionalPhase::Gas);
    const auto at350 = fluid.eos.mixingParameters(
        5.0e6, 350.0, MPMC::CompositionalPhase::Gas);
    const auto recover = [](const auto &m) {
        return 1.0 - m.Aij[0][1] / std::sqrt(m.Aij[0][0] * m.Aij[1][1]);
    };
    near(recover(at300), 0.0200, 2.0e-13, "mixing matrix BIP at reference T");
    near(recover(at350), 0.0300, 2.0e-13, "mixing matrix BIP away from reference T");

    auto callbackFluid = MPMC::cases::makeFluidSystem<Indices, CallbackFactoryConfig>();
    near(callbackFluid.eos.binaryInteractionCoefficient(0, 1, 350.0), 0.135, 1.0e-15,
         "factory arbitrary BIP callback must take precedence over linear coefficients");
}

void checkVolumeTranslationPreservesEquilibrium()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, ExtendedFactoryConfig>();
    require(fluid.eos.usesVolumeTranslation(), "factory must enable component volume translation");

    Eos unshifted = fluid.eos;
    unshifted.configureVolumeTranslation({0.0, 0.0});

    constexpr double p = 3.0e6;
    constexpr double T = 340.0;
    const Composition x{0.4, 0.6};
    const auto shiftedThermo = fluid.eos.phaseResult(
        p, T, x, MPMC::CompositionalPhase::Gas, false);
    const auto unshiftedThermo = unshifted.phaseResult(
        p, T, x, MPMC::CompositionalPhase::Gas, false);

    near(shiftedThermo.compressibility, unshiftedThermo.compressibility, 1.0e-15,
         "volume translation must not change cubic Z");
    for (std::size_t c = 0; c < x.size(); ++c)
    {
        near(shiftedThermo.fugacityCoefficient[c], unshiftedThermo.fugacityCoefficient[c], 1.0e-15,
             "volume translation must not change fugacity coefficients");
        near(shiftedThermo.fugacity[c], unshiftedThermo.fugacity[c], 1.0e-15,
             "volume translation must not change fugacity");
    }

    constexpr double R = MPMC::units::gasConstant;
    const double cMix = 0.4e-6 + 0.6 * 3.0e-6;
    const double expected = 1.0 /
        (shiftedThermo.compressibility * R * T / p - cMix);
    const double shiftedDensity = fluid.eos.molarDensity(
        p, T, x, shiftedThermo.compressibility);
    const double legacyDensity = unshifted.molarDensity(
        p, T, x, unshiftedThermo.compressibility);
    near(shiftedDensity, expected, 2.0e-14,
         "translated molar density must follow v=ZRT/P-sum(x_i c_i)");
    require(shiftedDensity > legacyDensity,
            "positive Peneloux-style c_i must reduce molar volume and raise density");

    const auto [mass, massDensity, viscosity] = fluid.eosFlowProperties(
        p, x, shiftedThermo.compressibility, T);
    (void)mass;
    (void)viscosity;
    const double molarMass = 0.4 * ExtendedFactoryConfig::Fluid::molarMass[0]
        + 0.6 * ExtendedFactoryConfig::Fluid::molarMass[1];
    near(massDensity, shiftedDensity * molarMass, 2.0e-13,
         "FluidSystem flow density must use translated EOS molar density");
}

void checkSwAqueousVolumeTranslationIsPhaseSpecific()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, ExtendedFactoryConfig>();
    Eos uncorrected = fluid.eos;

    Eos::SoreideWhitsonOptions base;
    base.waterComponent = 0;
    base.salinityMolality = 1.0;
    base.aqueousWaterBip[1] = [](double, double) { return 0.15; };
    uncorrected.configureSoreideWhitson(base);

    Eos corrected = uncorrected;
    base.aqueousVolumeTranslation = {0.0, 2.0e-6};
    corrected.configureSoreideWhitson(base);

    constexpr double p = 4.0e6;
    constexpr double T = 340.0;
    const Composition x{0.97, 0.03};
    const auto waterThermo = corrected.phaseResult(
        p, T, x, MPMC::CompositionalPhase::Water, false);
    const auto referenceThermo = uncorrected.phaseResult(
        p, T, x, MPMC::CompositionalPhase::Water, false);
    near(waterThermo.compressibility, referenceThermo.compressibility, 1.0e-15,
         "SW aqueous density correction must not change Z");
    for (std::size_t c = 0; c < x.size(); ++c)
    {
        near(waterThermo.fugacityCoefficient[c],
             referenceThermo.fugacityCoefficient[c], 1.0e-15,
             "SW aqueous density correction must not change fugacity coefficient");
    }

    const double rhoWater = corrected.molarDensity(
        p, T, x, waterThermo.compressibility, MPMC::CompositionalPhase::Water);
    const double rhoOil = corrected.molarDensity(
        p, T, x, waterThermo.compressibility, MPMC::CompositionalPhase::Oil);
    const double expectedWater = 1.0 /
        (1.0 / rhoOil - x[1] * 2.0e-6);
    near(rhoWater, expectedWater, 2.0e-14,
         "SW aqueous translation must apply only to the Water role");
    require(rhoWater > rhoOil,
            "positive aqueous CO2 translation must raise only water-phase density");
}



void checkVolumeTranslationAdDerivative()
{
    constexpr std::array<double, 2> tc{305.32, 469.7};
    constexpr std::array<double, 2> pc{4.872e6, 3.37e6};
    constexpr std::array<double, 2> vc{1.45e-4, 3.13e-4};
    constexpr std::array<double, 2> omega{0.0995, 0.251};
    constexpr std::array<double, 2> mw{0.03007, 0.07215};
    constexpr std::array<std::array<double, 2>, 2> kij{{{{0.0, 0.02}}, {{0.02, 0.0}}}};

    AdEos ad(
        0.45724, 0.07780,
        MPMC::CompositionalMixture<AdIndices>(tc, pc, vc, omega, mw, kij),
        1, 2.414213562373095, -0.414213562373095, 1.0e-30);
    ad.configureVolumeTranslation({1.0e-6, 3.0e-6});
    std::array<AdValue, 2> x{AdValue(0.4), AdValue(0.6)};
    const auto p = AdValue::createVariable(3.0e6, AdIndices::Primary::pressure);
    const AdValue T(340.0);
    const auto phase = ad.phaseResult(p, T, x, MPMC::CompositionalPhase::Gas, false);
    const auto rho = ad.molarDensity(p, T, x, phase.compressibility);
    const double adDerivative = rho.derivative(AdIndices::Primary::pressure);
    require(std::isfinite(adDerivative) && std::abs(adDerivative) > 0.0,
            "translated molar density must retain an AD pressure derivative");

    Eos scalar(
        0.45724, 0.07780,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        1, 2.414213562373095, -0.414213562373095, 1.0e-30);
    scalar.configureVolumeTranslation({1.0e-6, 3.0e-6});
    const Composition scalarX{0.4, 0.6};
    constexpr double dp = 50.0;
    const auto densityAt = [&](double pressure) {
        const auto z = scalar.phaseResult(
            pressure, 340.0, scalarX, MPMC::CompositionalPhase::Gas, false).compressibility;
        return scalar.molarDensity(pressure, 340.0, scalarX, z);
    };
    const double finiteDifference = (densityAt(3.0e6 + dp) - densityAt(3.0e6 - dp)) / (2.0 * dp);
    near(adDerivative, finiteDifference, 2.0e-5,
         "volume-translation AD pressure derivative must match finite difference");
}

void checkSwPrecedenceAndCpaTemperatureBip()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, ExtendedFactoryConfig>();
    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[1] = [](double, double) { return -0.075; };
    fluid.eos.configureSoreideWhitson(std::move(sw));
    fluid.eos.configureBinaryInteractionFunction(
        [](int i, int j, double T) {
            if (i == j)
                return 0.0;
            return 0.20 + 1.0e-3 * (T - 300.0);
        });

    constexpr double p = 4.0e6;
    constexpr double T = 340.0;
    const auto gas = fluid.eos.mixingParameters(p, T, MPMC::CompositionalPhase::Gas);
    const auto water = fluid.eos.mixingParameters(p, T, MPMC::CompositionalPhase::Water);
    const auto recover = [](const auto &m) {
        return 1.0 - m.Aij[0][1] / std::sqrt(m.Aij[0][0] * m.Aij[1][1]);
    };
    near(recover(gas), 0.24, 2.0e-13,
         "SW non-aqueous phase must honor generic temperature-dependent BIP");
    near(recover(water), -0.075, 2.0e-13,
         "SW aqueous water pair must retain aqueous-correlation precedence");

    // The same generic k_ij(T) hook also feeds the SRK cubic contribution of CPA.
    Eos cpa = fluid.eos;
    Eos::CubicPlusAssociationOptions options;
    constexpr double R = MPMC::units::gasConstant;
    for (int i = 0; i < 2; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        const double tc = ExtendedFactoryConfig::Fluid::criticalTemperature[idx];
        const double pc = ExtendedFactoryConfig::Fluid::criticalPressure[idx];
        const double omega = ExtendedFactoryConfig::Fluid::acentricFactor[idx];
        options.a0[idx] = 0.42748 * R * R * tc * tc / pc;
        options.b[idx] = 0.08664 * R * tc / pc;
        options.c1[idx] = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    }
    cpa.configureCubicPlusAssociation(options);
    cpa.configureBinaryInteractionFunction(
        [](int i, int j, double Tvalue) {
            return i == j ? 0.0 : 0.05 + 5.0e-4 * (Tvalue - 300.0);
        });
    Eos cpaConstant = cpa;
    cpaConstant.clearBinaryInteractionFunction();
    const Composition x{0.35, 0.65};
    const auto variable = cpa.phaseResult(p, T, x, MPMC::CompositionalPhase::Gas, false);
    const auto constant = cpaConstant.phaseResult(p, T, x, MPMC::CompositionalPhase::Gas, false);
    require(std::abs(variable.fugacityCoefficient[0] - constant.fugacityCoefficient[0]) > 1.0e-8,
            "CPA cubic contribution must consume the temperature-dependent BIP hook");
}

void checkChabab2019Co2AqueousBip()
{
    // Independent evaluation of Chabab et al. (2019), Eq. for the modified
    // aqueous CO2-H2O BIP. The fixed 304.13 K divisor is part of the fit.
    const auto reference = [](double temperature, double salinity) {
        const double tr = temperature / 304.13;
        return tr * (
            0.43575155 - 0.05766906744 * tr
            + 0.00826464849 * tr * salinity)
            + salinity * salinity * (0.00129539193 - 0.0016698848 * tr)
            - 0.47866096;
    };
    for (const auto &[temperature, salinity] : std::array{
             std::pair{303.15, 0.0},
             std::pair{323.15, 1.0},
             std::pair{373.15, 3.0},
             std::pair{423.15, 6.0}})
    {
        near(
            MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(
                temperature, salinity),
            reference(temperature, salinity), 2.0e-15,
            "Chabab 2019 aqueous CO2-H2O BIP identity");
    }

    bool rejected = false;
    try
    {
        (void)MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(
            323.15, -1.0);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    require(rejected, "Chabab 2019 BIP must reject negative salinity");
}

void checkArbitraryBipCallback()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, ExtendedFactoryConfig>();
    fluid.eos.configureBinaryInteractionFunction(
        [](int i, int j, double T) {
            if (i == j)
                return 0.0;
            return 0.15 + 1.0e-4 * (T - 320.0);
        });
    near(fluid.eos.binaryInteractionCoefficient(0, 1, 370.0), 0.155, 1.0e-15,
         "arbitrary cubic BIP callback");
}

void checkDefaultCrossAssociationBudget()
{
    // Public ThermoPack water/methanol sCPA parameters.  This liquid state is
    // a regression for the historical 100-iteration ceiling: the damped site
    // equations converge, but need more than 100 updates.
    Eos::CubicPlusAssociationOptions cpa;
    require(cpa.maximumAssociationIterations == 2000,
            "CPA default must retain the cross-association convergence budget");

    auto eos = makeWaterMethanolCpa<Indices>();
    eos.setThermodynamicProfilerEnabled(true);
    eos.resetThermodynamicProfile();
    const Composition x{0.185, 0.815};
    const auto liquid = eos.phaseResult(
        233183.46666657156, 389.375, x,
        MPMC::CompositionalPhase::Oil, false);
    require(std::isfinite(liquid.compressibility),
            "water/methanol cross-associating liquid state must converge by default");
    near(liquid.compressibility, 0.00292019189525, 2.0e-10,
         "water/methanol CPA liquid Z regression");
    const auto profile = eos.thermodynamicProfile();
    const double averageIterations =
        static_cast<double>(profile.cpaAssociationIterations) /
        static_cast<double>(profile.cpaAssociationIterativeCalls);
    require(averageIterations < 32.0,
            "general CPA Newton solve must not regress to slow fixed-point work");
}

void checkB2ReducedAssociationAdDerivative()
{
    auto ad = makeB2Cpa<AdIndices>();
    std::array<AdValue, 2> adComposition{AdValue(0.75), AdValue(0.25)};
    const auto pressure = AdValue::createVariable(
        8.0e6, AdIndices::Primary::pressure);
    const auto phase = ad.phaseResult(
        pressure, AdValue(305.0), adComposition,
        MPMC::CompositionalPhase::Water, false);

    const double zDerivative =
        phase.compressibility.derivative(AdIndices::Primary::pressure);
    const double phiDerivative =
        phase.fugacityCoefficient[0].derivative(AdIndices::Primary::pressure);
    require(std::isfinite(zDerivative) && std::isfinite(phiDerivative),
            "reduced CPA association solve must retain finite AD derivatives");

    auto scalar = makeB2Cpa<Indices>();
    const Composition composition{0.75, 0.25};
    constexpr double dp = 100.0;
    const auto evaluate = [&](double p) {
        return scalar.phaseResult(
            p, 305.0, composition,
            MPMC::CompositionalPhase::Water, false);
    };
    const auto lower = evaluate(8.0e6 - dp);
    const auto upper = evaluate(8.0e6 + dp);
    const double zFiniteDifference =
        (upper.compressibility - lower.compressibility) / (2.0 * dp);
    const double phiFiniteDifference =
        (upper.fugacityCoefficient[0] - lower.fugacityCoefficient[0]) /
        (2.0 * dp);
    const auto derivativeNear = [](double actual, double expected) {
        return std::abs(actual - expected) <=
            2.0e-4 * std::max(std::abs(expected), 1.0e-16);
    };
    require(derivativeNear(zDerivative, zFiniteDifference),
            "reduced CPA Z pressure derivative must match finite difference");
    require(derivativeNear(phiDerivative, phiFiniteDifference),
            "reduced CPA fugacity derivative must match finite difference");
}

void checkGeneralAssociationAdDerivative()
{
    auto ad = makeWaterMethanolCpa<AdIndices>();
    std::array<AdValue, 2> adComposition{AdValue(0.185), AdValue(0.815)};
    const auto pressure = AdValue::createVariable(
        233183.46666657156, AdIndices::Primary::pressure);
    const auto phase = ad.phaseResult(
        pressure, AdValue(389.375), adComposition,
        MPMC::CompositionalPhase::Oil, false);
    const double zDerivative =
        phase.compressibility.derivative(AdIndices::Primary::pressure);
    const double phiDerivative =
        phase.fugacityCoefficient[0].derivative(AdIndices::Primary::pressure);
    require(std::isfinite(zDerivative) && std::isfinite(phiDerivative),
            "general CPA Newton solve must retain finite AD derivatives");

    auto scalar = makeWaterMethanolCpa<Indices>();
    const Composition composition{0.185, 0.815};
    constexpr double dp = 10.0;
    const auto evaluate = [&](double p) {
        return scalar.phaseResult(
            p, 389.375, composition,
            MPMC::CompositionalPhase::Oil, false);
    };
    const auto lower = evaluate(233183.46666657156 - dp);
    const auto upper = evaluate(233183.46666657156 + dp);
    const double zFiniteDifference =
        (upper.compressibility - lower.compressibility) / (2.0 * dp);
    const double phiFiniteDifference =
        (upper.fugacityCoefficient[0] - lower.fugacityCoefficient[0]) /
        (2.0 * dp);
    const auto derivativeNear = [](double actual, double expected) {
        return std::abs(actual - expected) <=
            5.0e-4 * std::max(std::abs(expected), 1.0e-16);
    };
    require(derivativeNear(zDerivative, zFiniteDifference),
            "general CPA Z pressure derivative must match finite difference");
    require(derivativeNear(phiDerivative, phiFiniteDifference),
            "general CPA fugacity derivative must match finite difference");

    const int compositionIndex = AdIndices::Primary::liquidComposition[0];
    const auto waterFraction = AdValue::createVariable(0.185, compositionIndex);
    std::array<AdValue, 2> variableComposition{
        waterFraction, AdValue(1.0) - waterFraction};
    const auto compositionPhase = ad.phaseResult(
        AdValue(233183.46666657156), AdValue(389.375), variableComposition,
        MPMC::CompositionalPhase::Oil, false);
    const double zCompositionDerivative =
        compositionPhase.compressibility.derivative(compositionIndex);
    const double phiCompositionDerivative =
        compositionPhase.fugacityCoefficient[0].derivative(compositionIndex);
    constexpr double dx = 1.0e-6;
    const auto compositionAt = [&](double xWater) {
        return scalar.phaseResult(
            233183.46666657156, 389.375,
            Composition{xWater, 1.0 - xWater},
            MPMC::CompositionalPhase::Oil, false);
    };
    const auto compositionLower = compositionAt(0.185 - dx);
    const auto compositionUpper = compositionAt(0.185 + dx);
    const double zCompositionFiniteDifference =
        (compositionUpper.compressibility - compositionLower.compressibility) /
        (2.0 * dx);
    const double phiCompositionFiniteDifference =
        (compositionUpper.fugacityCoefficient[0] -
         compositionLower.fugacityCoefficient[0]) / (2.0 * dx);
    require(derivativeNear(
                zCompositionDerivative, zCompositionFiniteDifference),
            "general CPA Z composition derivative must match finite difference");
    require(derivativeNear(
                phiCompositionDerivative, phiCompositionFiniteDifference),
            "general CPA fugacity composition derivative must match finite difference");
}

void checkIapwsGarciaAqueousVolume()
{
    constexpr double pressure = 5.16e6;
    constexpr double temperature = 333.15;
    const double pureWaterDensity =
        MPMC::IapwsIf97Region1::density(pressure, temperature);
    near(pureWaterDensity, 985.4040020947351, 2.0e-12,
         "IAPWS-IF97 Region-1 main-state water density");
    near(MPMC::IapwsIf97Region2::density(30.0e6, 700.0),
         184.18016875974072, 2.0e-12,
         "IAPWS-IF97 Table-15 Region-2 water density");
    near(MPMC::IapwsIf97WaterDensity::density(24.0e6, 673.15),
         148.56108348732423, 2.0e-12,
         "IAPWS-IF97 supercritical-water benchmark density");
    near(MPMC::IapwsIf97Region3::pressureFromDensity(500.0, 650.0),
         25.5837018e6, 2.0e-9,
         "IAPWS-IF97 Table-33 Region-3 pressure");
    near(MPMC::IapwsIf97Region3::density(25.5837018e6, 650.0),
         500.0, 5.0e-9,
         "IAPWS-IF97 Region-3 density inversion");
    const double co2Volume =
        MPMC::Garcia2001Co2ApparentVolume::molarVolume(temperature);
    near(co2Volume, 34.7964496e-6, 2.0e-14,
         "Garcia aqueous CO2 apparent molar volume at 60 C");

    MPMC::IapwsGarciaAqueousVolume<2> model(0, 1, 0.018015268);
    const Composition x{0.9873, 0.0127};
    const double molarVolume = model.molarVolume(pressure, temperature, x);
    const double mixtureMolarMass =
        x[0] * 0.018015268 + x[1] * 0.0440098;
    near(mixtureMolarMass / molarVolume, 992.080, 2.0e-6,
         "IAPWS-Garcia dissolved-CO2 aqueous density regression");

    MPMC::IapwsGarciaAqueousVolume<2> adModel(0, 1, 0.018015268);
    const auto adPressure = AdValue::createVariable(
        pressure, AdIndices::Primary::pressure);
    std::array<AdValue, 2> adComposition{AdValue(x[0]), AdValue(x[1])};
    const auto adVolume = adModel.molarVolume(
        adPressure, AdValue(temperature), adComposition);
    const double adDerivative = adVolume.derivative(AdIndices::Primary::pressure);
    constexpr double dp = 100.0;
    const double finiteDifference =
        (model.molarVolume(pressure + dp, temperature, x) -
         model.molarVolume(pressure - dp, temperature, x)) / (2.0 * dp);
    near(adDerivative, finiteDifference, 2.0e-6,
         "IAPWS-Garcia pressure derivative must match finite difference");

    const auto supercriticalPressure = AdValue::createVariable(
        24.0e6, AdIndices::Primary::pressure);
    const auto supercriticalDensity = MPMC::IapwsIf97WaterDensity::density(
        supercriticalPressure, AdValue(673.15));
    constexpr double supercriticalDp = 100.0;
    const double supercriticalFiniteDifference =
        (MPMC::IapwsIf97WaterDensity::density(
             24.0e6 + supercriticalDp, 673.15) -
         MPMC::IapwsIf97WaterDensity::density(
             24.0e6 - supercriticalDp, 673.15)) /
        (2.0 * supercriticalDp);
    near(supercriticalDensity.derivative(AdIndices::Primary::pressure),
         supercriticalFiniteDifference, 2.0e-6,
         "IAPWS Region-2 pressure derivative must match finite difference");

    const auto densePressure = AdValue::createVariable(
        36.0e6, AdIndices::Primary::pressure);
    const auto denseDensity = MPMC::IapwsIf97WaterDensity::density(
        densePressure, AdValue(673.15));
    constexpr double denseDp = 100.0;
    const double denseFiniteDifference =
        (MPMC::IapwsIf97WaterDensity::density(36.0e6 + denseDp, 673.15) -
         MPMC::IapwsIf97WaterDensity::density(36.0e6 - denseDp, 673.15)) /
        (2.0 * denseDp);
    near(denseDensity.derivative(AdIndices::Primary::pressure),
         denseFiniteDifference, 2.0e-6,
         "IAPWS Region-3 pressure derivative must match finite difference");

    MPMC::IapwsGarciaAqueousVolume<3> boundedModel(
        0, 1, 0.018015268, 1.0e-4);
    require(boundedModel.compositionSupported({0.98, 0.01995, 5.0e-5}),
            "IAPWS-Garcia domain must admit trace third solute");
    require(!boundedModel.compositionSupported({0.80, 0.02, 0.18}),
            "IAPWS-Garcia domain must reject a hydrocarbon-rich Water-role basin");
}

void checkMcBrideWrightAqueousViscosity()
{
    MPMC::McBrideWright2015AqueousViscosity<2> model(0, 1);
    constexpr double pressure = 5.16e6;
    constexpr double temperature = 333.15;
    near(model.viscosity(pressure, temperature, Composition{1.0, 0.0}),
         0.469091e-3, 3.0e-6,
         "McBride-Wright pure-water main-state viscosity");
    near(model.viscosity(
             pressure, temperature, Composition{1.0 - 0.01304222, 0.01304222}),
         0.483480e-3, 3.0e-6,
         "McBride-Wright dissolved-CO2 viscosity regression");
}

void checkIapws2008AqueousViscosity()
{
    MPMC::Iapws2008IndustrialAqueousViscosity<2> model(0, 0.02);
    near(model.viscosity(
             0.101325e6, 298.15, Composition{1.0, 0.0}),
         889.735100e-6, 2.0e-6,
         "IAPWS-2008 ambient-water verification point");
    const double supercritical = model.viscosity(
        27.74e6, 653.2, Composition{0.9966, 0.0034});
    require(std::isfinite(supercritical) && supercritical > 0.0 &&
                supercritical < 2.0e-4,
            "IAPWS-2008 supercritical-water viscosity must be finite and plausible");
    bool rejected = false;
    try
    {
        (void)model.viscosity(
            27.74e6, 653.2, Composition{0.95, 0.05});
    }
    catch (const std::runtime_error &)
    {
        rejected = true;
    }
    require(rejected,
            "IAPWS-2008 pure-water closure must reject an unsupported solute fraction");
}

void checkAqueousCompositionDomain()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, ExtendedFactoryConfig>();
    fluid.eos.configureAqueousCompositionDomain(0, 0.02);
    require(fluid.eos.aqueousVolumeCompositionSupported({0.99, 0.01}),
            "validated water-rich composition must remain in the aqueous-property domain");
    require(!fluid.eos.aqueousVolumeCompositionSupported({0.97, 0.03}),
            "hydrocarbon-rich composition must not select aqueous-only properties");
}

void checkPrCpaDryLimit()
{
    constexpr std::array<double, 2> tc{304.1282, 617.70};
    constexpr std::array<double, 2> pc{7.37703e6, 2.110e6};
    constexpr std::array<double, 2> vc{9.41185e-5, 6.03e-4};
    constexpr std::array<double, 2> omega{0.2249, 0.492328};
    constexpr std::array<double, 2> mw{0.0440095, 0.1422817};
    constexpr std::array<std::array<double, 2>, 2> kij{{
        {{0.0, 0.1141}}, {{0.1141, 0.0}}}};
    Eos pr(
        0.45724, 0.07780,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        5, 2.414213562373095, -0.414213562373095, 1.0e-30);
    Eos prCpa = pr;
    Eos::CubicPlusAssociationOptions options;
    constexpr double R = MPMC::units::gasConstant;
    for (std::size_t c = 0; c < 2; ++c)
    {
        options.a0[c] = 0.45724 * R * R * tc[c] * tc[c] / pc[c];
        options.b[c] = 0.07780 * R * tc[c] / pc[c];
        options.c1[c] =
            0.37464 + 1.54226 * omega[c] - 0.26992 * omega[c] * omega[c];
        if (omega[c] > 0.49)
        {
            options.c1[c] = 0.379642 + 1.48503 * omega[c]
                - 0.164423 * omega[c] * omega[c]
                + 0.016666 * omega[c] * omega[c] * omega[c];
        }
    }
    options.physicalTerm = MPMC::CpaCubicPhysicalTerm::PengRobinson;
    prCpa.configureCubicPlusAssociation(options);
    const Composition x{0.4166688, 1.0 - 0.4166688};
    for (const bool liquid : {false, true})
    {
        const auto role = liquid ? MPMC::CompositionalPhase::Oil
                                 : MPMC::CompositionalPhase::Gas;
        const auto reference = pr.phaseResult(5.16e6, 333.15, x, role, false);
        const auto actual = prCpa.phaseResult(5.16e6, 333.15, x, role, false);
        near(actual.compressibility, reference.compressibility, 2.0e-8,
             "non-associating PR-CPA must recover PR Z");
        for (std::size_t c = 0; c < x.size(); ++c)
            near(actual.fugacityCoefficient[c], reference.fugacityCoefficient[c],
                 3.0e-8, "non-associating PR-CPA must recover PR fugacity");
    }
}

} // namespace

int main()
{
    checkFactoryTemperatureDependentBip();
    checkVolumeTranslationPreservesEquilibrium();
    checkSwAqueousVolumeTranslationIsPhaseSpecific();
    checkArbitraryBipCallback();
    checkSwPrecedenceAndCpaTemperatureBip();
    checkChabab2019Co2AqueousBip();
    checkVolumeTranslationAdDerivative();
    checkDefaultCrossAssociationBudget();
    checkB2ReducedAssociationAdDerivative();
    checkGeneralAssociationAdDerivative();
    checkIapwsGarciaAqueousVolume();
    checkMcBrideWrightAqueousViscosity();
    checkIapws2008AqueousViscosity();
    checkAqueousCompositionDomain();
    checkPrCpaDryLimit();
    std::cout << "thermo_empirical_parameters_test: PASS\n";
    return 0;
}
