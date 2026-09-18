/**
 * @file thermodynamic_failure_modes_test.cpp
 * @brief 回归：覆盖 PR/SW/CPA flash 的失败恢复、稳定性证书与 BIP 契约。
 */
#include "../../../case/three_eos_3d_compare/case_config.hpp"

#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;
using Composition = std::array<double, Indices::numComponents>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void checkMaterialClosure(const Composition &z, const Flash::Result &result)
{
    double betaSum = 0.0;
    Composition reconstructed{};
    for (int p = 0; p < 3; ++p)
    {
        const std::size_t phase = static_cast<std::size_t>(p);
        betaSum += result.phaseMoleFraction[phase];
        for (int c = 0; c < Indices::numComponents; ++c)
        {
            reconstructed[static_cast<std::size_t>(c)] +=
                result.phaseMoleFraction[phase] *
                result.composition[phase][static_cast<std::size_t>(c)];
        }
    }
    require(std::abs(betaSum - 1.0) < 2.0e-8,
            "flash phase fractions must close");
    for (int c = 0; c < Indices::numComponents; ++c)
    {
        require(std::abs(reconstructed[static_cast<std::size_t>(c)] -
                         z[static_cast<std::size_t>(c)]) < 2.0e-8,
                "flash component material balance must close");
    }
}

template <class FactoryConfig>
Flash makeFlash()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    return Flash(fluid.eos, options);
}

template <class FactoryConfig>
void requireStableFlash(
    double pressure,
    const Composition &z,
    const std::string &label)
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    Flash flash(fluid.eos, options);

    Composition normalized = z;
    double sum = 0.0;
    for (double value : normalized)
        sum += value;
    for (double &value : normalized)
        value /= sum;

    const auto result = flash.flash(
        pressure, CaseConfig::InitialState::temperature, normalized);
    require(result.converged, label + " unrestricted flash must converge");
    checkMaterialClosure(normalized, result);

    if (result.presence.count() < 3)
    {
        const auto stability = flash.stabilityTest(
            pressure,
            CaseConfig::InitialState::temperature,
            normalized,
            result.presence,
            result.composition);
        require(stability.valid,
                label + " reduced-set result must carry a valid stability certificate");
        require(stability.stable,
                label + " reduced-set result must be thermodynamically stable");
    }
}

void checkPreviouslyMissedReducedSets()
{
    // v67.1 default PR path failed at this P-T-z although O+W is a stable
    // restricted equilibrium.  v67.2 must find that branch automatically.
    const Composition prState{
        0.05009449, 0.16662213, 0.33957218, 0.40901014, 0.03470106};
    requireStableFlash<CaseConfig::PrFactoryConfig>(
        162.9742365 * CaseConfig::bar,
        prState,
        "PR reduced-set recovery");

    // A deterministic CPA state from the v67.1 random robustness audit.  The
    // direct unrestricted path failed, while a stable G+W branch existed.
    const Composition cpaState{
        0.0639856, 0.0929489, 0.171507, 0.267589, 0.403969};
    requireStableFlash<CaseConfig::CpaFactoryConfig>(
        64.5192 * CaseConfig::bar,
        cpaState,
        "CPA reduced-set recovery");
}

template <class FactoryConfig>
void deterministicRobustnessSweep(const std::string &label)
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    Flash flash(fluid.eos, options);

    std::mt19937_64 generator(1234567);
    std::uniform_real_distribution<double> compositionSample(0.001, 1.0);
    std::uniform_real_distribution<double> pressureBar(20.0, 180.0);

    // This sweep is intentionally small enough for the default unit suite but
    // traverses the broad simplex that exposed v67.1 path-dependent failures.
    constexpr int samples = 48;
    for (int sample = 0; sample < samples; ++sample)
    {
        Composition z{};
        double sum = 0.0;
        for (double &value : z)
        {
            value = compositionSample(generator);
            sum += value;
        }
        for (double &value : z)
            value /= sum;

        const double pressure = pressureBar(generator) * CaseConfig::bar;
        const auto result = flash.flash(
            pressure, CaseConfig::InitialState::temperature, z);
        require(result.converged,
                label + " deterministic robustness sweep contains a failed flash");
        checkMaterialClosure(z, result);

        if (result.presence.count() < 3)
        {
            const auto stability = flash.stabilityTest(
                pressure,
                CaseConfig::InitialState::temperature,
                z,
                result.presence,
                result.composition);
            require(stability.valid && stability.stable,
                    label + " robustness sweep returned an uncertified reduced phase set");
        }
    }
}


void checkKerogenCpaCoincidentNonaqueousOnsetRecovery()
{
    using LocalMixture = MPMC::CompositionalMixture<Indices>;
    using LocalEos = MPMC::CubicEquationOfState<Indices>;

    constexpr std::array<double, 5> tc{
        647.096, 515.231, 736.423, 870.437, 982.864};
    constexpr std::array<double, 5> pc{
        22.064e6, 3.192179e6, 1.751155e6, 1.138075e6, 0.808714e6};
    constexpr std::array<double, 5> vc{
        55.948074534e-6, 358.490e-6, 824.073e-6,
        1318.853e-6, 1793.890e-6};
    constexpr std::array<double, 5> omega{
        0.3443, 0.269328, 0.583473, 0.896554, 1.215676};
    constexpr std::array<double, 5> mw{
        0.01801528, 0.082392, 0.215583, 0.387349, 0.660132};
    constexpr std::array<std::array<double, 5>, 5> kij{{
        {{0.0, 0.044, -0.054, 0.0, 0.0}},
        {{0.044, 0.0, 0.0, 0.0, 0.0}},
        {{-0.054, 0.0, 0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0, 0.0, 0.0}}
    }};

    LocalEos eos(
        0.42748, 0.08664,
        LocalMixture(tc, pc, vc, omega, mw, kij),
        1, 1.0, 0.0, 1.0e-30);

    LocalEos::CubicPlusAssociationOptions cpa;
    cpa.a0 = {
        0.12277, 2.45754175486, 9.15196330734,
        19.6737874211, 35.3001105724};
    cpa.b = {
        1.4515e-5, 1.16269921130e-4, 3.02939137177e-4,
        5.50958755266e-4, 8.75489809612e-4};
    cpa.c1 = {
        0.67359, 0.891155659401, 1.33846893146,
        1.74970559881, 2.13336923189};
    cpa.associationEnergy[0] = 16655.0;
    cpa.associationVolume[0] = 0.0692;
    cpa.donorSites[0] = 2;
    cpa.acceptorSites[0] = 2;
    cpa.physicalTerm = MPMC::CpaCubicPhysicalTerm::SoaveRedlichKwong;
    cpa.radialDistribution = MPMC::CpaRadialDistribution::Simplified;
    eos.configureCubicPlusAssociation(cpa);

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    options.maximumIterations = 180;
    options.maximumStabilityIterations = 120;
    Flash flash(eos, options);

    const Composition z{
        0.7610625,
        0.008680599375,
        0.097192606875,
        0.0777552801875,
        0.0553090135625};

    for (double pressure : {26.0e6, 26.25e6})
    {
        const auto result = flash.flash(pressure, 653.15, z);
        require(result.converged,
                "kerogen CPA onset regression must converge");
        require(result.presence.contains(MPMC::CompositionalPhase::Oil) &&
                    result.presence.contains(MPMC::CompositionalPhase::Water) &&
                    !result.presence.contains(MPMC::CompositionalPhase::Gas),
                "kerogen CPA onset regression must recover stable O+W");
        checkMaterialClosure(z, result);

        const auto stability = flash.stabilityTest(
            pressure, 653.15, z, result.presence, result.composition);
        require(stability.valid && stability.stable,
                "kerogen CPA O+W onset state must carry a stable certificate");
        require(result.composition[2][0] > result.composition[0][0],
                "kerogen CPA Water role must remain more H2O-rich than Oil");
    }
}


void checkStaticBipSymmetryContract()
{
    using Mixture = MPMC::CompositionalMixture<Indices>;
    Mixture::ScalarArray tc{};
    Mixture::ScalarArray pc{};
    Mixture::ScalarArray vc{};
    Mixture::ScalarArray omega{};
    Mixture::ScalarArray molarMass{};
    Mixture::BinaryInteractionMatrix bip{};

    for (int i = 0; i < Indices::numComponents; ++i)
    {
        const auto c = static_cast<std::size_t>(i);
        tc[c] = 300.0 + 10.0 * i;
        pc[c] = 4.0e6 + 1.0e5 * i;
        vc[c] = 1.0e-4 + 1.0e-5 * i;
        omega[c] = 0.1 + 0.01 * i;
        molarMass[c] = 0.02 + 0.01 * i;
    }
    bip[0][1] = 0.1;
    bip[1][0] = 0.2;

    bool rejected = false;
    try
    {
        (void)Mixture(tc, pc, vc, omega, molarMass, bip);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    require(rejected,
            "static BIP matrix must reject asymmetric k_ij before EOS evaluation");
}

void checkTemperatureDependentBipSymmetryContract()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::PrFactoryConfig>();
    fluid.eos.configureBinaryInteractionFunction(
        [](int i, int j, double) {
            if (i == j)
                return 0.0;
            return i < j ? 0.1 : 0.2;
        });

    bool rejected = false;
    try
    {
        (void)fluid.eos.mixingParameters(
            CaseConfig::InitialState::pressure,
            CaseConfig::InitialState::temperature,
            MPMC::CompositionalPhase::Oil);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    require(rejected,
            "temperature-dependent BIP callback must reject asymmetric k_ij(T)");
}

} // namespace

int main()
{
    try
    {
        checkPreviouslyMissedReducedSets();
        checkKerogenCpaCoincidentNonaqueousOnsetRecovery();
        deterministicRobustnessSweep<CaseConfig::PrFactoryConfig>("PR");
        deterministicRobustnessSweep<CaseConfig::SwFactoryConfig>("SW");
        deterministicRobustnessSweep<CaseConfig::CpaFactoryConfig>("CPA");
        checkStaticBipSymmetryContract();
        checkTemperatureDependentBipSymmetryContract();
        std::cout << "Thermodynamic failure-mode regression: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Thermodynamic failure-mode regression: FAIL: "
                  << e.what() << '\n';
        return 1;
    }
}
