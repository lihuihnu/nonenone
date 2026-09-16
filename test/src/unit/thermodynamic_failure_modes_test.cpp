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
