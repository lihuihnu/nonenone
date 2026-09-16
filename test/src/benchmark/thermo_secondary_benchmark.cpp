/**
 * @file thermo_secondary_benchmark.cpp
 * @brief 性能基准：测量 `thermo_secondary_benchmark` 对应核心路径的计算开销。
 */
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/state/three_phase_equilibrium.hpp>

#include <array>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <utility>

namespace
{
using Config = MPMC::CompositionalModelConfig<
    4,
    true,
    false,
    false,
    false,
    false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;

Eos makeEos()
{
    MPMC::CompositionalMixture<Indices> mixture(
        {304.2, 190.6, 717.0, 647.3},
        {73.8e5, 46.0e5, 14.2e5, 220.5e5},
        {9.4e-5, 9.9e-5, 9.0e-4, 5.6e-5},
        {0.225, 0.008, 0.742, 0.344},
        {0.0440098, 0.016043, 0.22644, 0.01801528},
        {{0.0, 0.1000, 0.1250, 0.1896},
         {0.1000, 0.0, 0.0780, 0.4850},
         {0.1250, 0.0780, 0.0, 0.5000},
         {0.1896, 0.4850, 0.5000, 0.0}});

    return Eos(
        0.4572355,
        0.0779691,
        std::move(mixture),
        1,
        2.414213562373095,
        -0.414213562373095,
        1.0e-30);
}

MPMC::FluidSystem<Indices> makeFluid()
{
    MPMC::FluidSystem<Indices> fluid(
        {800.0, 20.0, 1000.0},
        {1.0e-3, 1.0e-5, 1.0e-4},
        makeEos(),
        {"CO2", "CH4", "nC16", "H2O"},
        {"Oil", "Gas", "Water"},
        350.0);
    fluid.configureFullyCompositionalThreePhase(3);
    return fluid;
}
} // namespace

int main()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    constexpr double pressure = 50.0e5;
    constexpr double temperature = 350.0;
    const std::array<double, Indices::numComponents> overall{
        0.75, 0.025, 0.025, 0.20};

    const auto flash = equilibrium.flashPTZ(pressure, temperature, overall);
    if (!flash.converged)
        return 2;

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = pressure;
    MPMC::PhaseStateData<Indices> phaseState{};
    equilibrium.assignFlashResult(primary, phaseState, flash);

    constexpr std::size_t iterations = 200000;
    double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < iterations; ++i)
    {
        equilibrium.updateSecondary(primary, phaseState);
        checksum += phaseState.phaseCompressibility[0];
    }
    const auto stop = std::chrono::steady_clock::now();

    const double elapsedMs =
        std::chrono::duration<double, std::milli>(stop - start).count();
    const double nsPerUpdate = elapsedMs * 1.0e6 /
        static_cast<double>(iterations);

    std::cout << std::fixed << std::setprecision(3)
              << "iterations=" << iterations
              << " elapsed_ms=" << elapsedMs
              << " ns_per_update=" << nsPerUpdate
              << " checksum=" << std::setprecision(12) << checksum
              << '\n';
    return 0;
}
