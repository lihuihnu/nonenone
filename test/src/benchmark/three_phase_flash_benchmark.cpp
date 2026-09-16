/**
 * @file three_phase_flash_benchmark.cpp
 * @brief 性能基准：测量 `three_phase_flash_benchmark` 对应核心路径的计算开销。
 */
#include "../../../case/3p4c_pr_reservoir/case_config.hpp"

#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/state/three_phase_equilibrium.hpp>

#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    false,
    false,
    false,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;
} // namespace

int main()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    constexpr std::size_t iterations = 500;
    double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < iterations; ++i)
    {
        const auto result = equilibrium.flashPTZ(
            CaseConfig::InitialState::pressure,
            CaseConfig::InitialState::temperature,
            CaseConfig::InitialState::overallComposition);
        if (!result.converged)
            return 2;
        checksum += result.phaseMoleFraction[0] +
            result.saturation[2] +
            result.composition[1][0];
    }
    const auto stop = std::chrono::steady_clock::now();

    const double elapsedMs =
        std::chrono::duration<double, std::milli>(stop - start).count();
    const double microsecondsPerFlash =
        elapsedMs * 1000.0 / static_cast<double>(iterations);

    std::cout << std::fixed << std::setprecision(3)
              << "iterations=" << iterations
              << " elapsed_ms=" << elapsedMs
              << " us_per_flash=" << microsecondsPerFlash
              << " checksum=" << std::setprecision(15) << checksum
              << '\n';
    return 0;
}
