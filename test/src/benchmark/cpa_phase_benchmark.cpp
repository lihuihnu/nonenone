/**
 * @file cpa_phase_benchmark.cpp
 * @brief 性能基准：测量 `cpa_phase_benchmark` 对应核心路径的计算开销。
 */
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/cubic_eos.hpp>

#include <array>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <utility>

namespace
{
using Config = MPMC::CompositionalModelConfig<
    2,
    true,
    false,
    false,
    false,
    false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;

Eos makeCpa()
{
    MPMC::CompositionalMixture<Indices> mixture(
        {647.096, 190.564},
        {220.64e5, 45.99e5},
        {5.6e-5, 9.9e-5},
        {0.344, 0.011},
        {0.01801528, 0.016043},
        {{{0.0, 0.0}, {0.0, 0.0}}});

    Eos eos(
        0.45724,
        0.07780,
        std::move(mixture),
        1,
        2.414213562373095,
        -0.414213562373095);

    Eos::CubicPlusAssociationOptions cpa;
    cpa.a0[0] = 0.12277;
    cpa.b[0] = MPMC::StandardCpaWater4C::b;
    cpa.c1[0] = 0.6736;
    cpa.associationEnergy[0] = 16655.0;
    cpa.associationVolume[0] = 0.0692;
    cpa.donorSites[0] = 2;
    cpa.acceptorSites[0] = 2;

    constexpr double R = MPMC::units::gasConstant;
    constexpr double tc = 190.564;
    constexpr double pc = 45.99e5;
    constexpr double omega = 0.011;
    cpa.a0[1] = 0.42748 * R * R * tc * tc / pc;
    cpa.b[1] = 0.08664 * R * tc / pc;
    cpa.c1[1] = 0.480 + 1.574 * omega - 0.176 * omega * omega;
    eos.configureCubicPlusAssociation(std::move(cpa));
    eos.configureCpaTemperatureCache(300.0);
    return eos;
}

Eos makeWaterMethanolCpa()
{
    MPMC::CompositionalMixture<Indices> mixture(
        {647.3, 512.6},
        {22.0483e6, 8.0959e6},
        {5.60e-5, 1.18e-4},
        {0.344, 0.559},
        {0.018015, 0.032042},
        {{{0.0, -0.09}, {-0.09, 0.0}}});
    Eos eos(0.42748, 0.08664, std::move(mixture), 1, 1.0, 0.0, 1.0e-30);
    Eos::CubicPlusAssociationOptions cpa;
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

void runBenchmark(
    std::string_view name,
    Eos &eos,
    double pressure,
    double temperature,
    const std::array<double, Indices::numComponents> &composition,
    MPMC::CompositionalPhase phase,
    std::size_t iterations)
{
    eos.setThermodynamicProfilerEnabled(true);
    eos.resetThermodynamicProfile();
    double checksum = 0.0;
    double lastCompressibility = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < iterations; ++i)
    {
        const auto result = eos.phaseResult(
            pressure, temperature, composition, phase, false);
        lastCompressibility = result.compressibility;
        checksum += result.compressibility + result.fugacityCoefficient[0];
    }
    const auto stop = std::chrono::steady_clock::now();
    const auto profile = eos.thermodynamicProfile();
    const double elapsedMs =
        std::chrono::duration<double, std::milli>(stop - start).count();
    const double iterationsPerAssociation = profile.cpaAssociationCalls > 0
        ? static_cast<double>(profile.cpaAssociationIterations) /
            static_cast<double>(profile.cpaAssociationCalls)
        : 0.0;
    std::cout << std::fixed << std::setprecision(3)
              << "name=" << name
              << " phases=" << iterations
              << " elapsed_ms=" << elapsedMs
              << " us_per_phase="
              << elapsedMs * 1000.0 / static_cast<double>(iterations)
              << " site_iterations_per_call=" << iterationsPerAssociation
              << " association_calls=" << profile.cpaAssociationCalls
              << " compressibility=" << std::setprecision(12)
              << lastCompressibility
              << " checksum=" << std::setprecision(12) << checksum
              << '\n';
}
} // namespace

int main()
{
    auto water = makeCpa();
    runBenchmark(
        "water_only", water, 100.0e5, 300.0, {0.95, 0.05},
        MPMC::CompositionalPhase::Water, 5000);

    auto waterMethanol = makeWaterMethanolCpa();
    runBenchmark(
        "water_methanol", waterMethanol, 233183.46666657156, 389.375,
        {0.185, 0.815}, MPMC::CompositionalPhase::Oil, 2000);
    return 0;
}
